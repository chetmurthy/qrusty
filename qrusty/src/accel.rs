use num_complex::Complex64;
use std::time::Instant;

#[derive(Debug, PartialEq)]
pub struct UnsafeVectors {
    pub data : Vec<Complex64>,
    pub indices: Vec<u64>,
    pub indptr: Vec<u64>,
}

pub fn make_unsafe_vectors(z: &Vec<bool>,
                      x: &Vec<bool>,
                      coeff: Complex64,
		      phase: i64,
		      group_phase: bool,
                ) -> std::result::Result<UnsafeVectors, &'static str> {
    let debug = false ;
    let timings = false ;
    let now = Instant::now();

    if z.len() != x.len() {
       Err("z and x have differing lengths")
    }
    else {
        let num_qubits = z.len() ;
	let mut mut_phase = phase ;

	if debug { println!("1: z={:?} x={:?} num_qubits={} mut_phase={}", z, x, num_qubits, mut_phase) ; }

	if group_phase {
	    let dotprod : i64 = x.iter()
                .zip(z.iter())
                .filter(|(x,z)| **x && **z)
                .count() as i64 ;
	    if debug { println!("2: dotprod={}", dotprod) ; }
	    mut_phase += dotprod ;
	    mut_phase = mut_phase % 4 ;
	}

	if debug { println!("2: mut_phase={}", mut_phase) ; }

	let dim =  1 << num_qubits ;
        let twos_array : Vec<u64> =
            (0..num_qubits)
            .map(|i| 1<<i)
            .collect() ;

	if debug { println!("3: twos_array={:?}", twos_array) ; }

        let x_indices : u64 =
            x.iter().enumerate()
            .filter(|(_,x)| **x)
            .map(|(i,_)| twos_array[i])
            .sum() ;

        let z_indices : u64 =
            z.iter().enumerate()
            .filter(|(_,z)| **z)
            .map(|(i,_)| twos_array[i])
            .sum() ;

	if debug { println!("4: x_indices={} z_indices={}", x_indices, z_indices) ; }


        if timings { println!("BEFORE indptr: {} ms", now.elapsed().as_millis()); }


        let vecsize = dim+1 ;
        let indptr : Vec<u64> = (0..vecsize).collect() ;

        if timings { println!("BEFORE indices: {} ms", now.elapsed().as_millis()); }

        let indices =
            (&indptr[..dim as usize]).iter()
            .map(|ind| ind ^ x_indices)
            .collect() ;

	let coeff = match phase % 4 {
	    0 => Complex64::new(1.0, 0.0) * coeff,
	    1 => Complex64::new(0.0, -1.0) * coeff,
	    2 => Complex64::new(-1.0, 0.0) * coeff,
	    3 => Complex64::new(0.0, 1.0) * coeff,
	    _ => coeff // really should be assert!(false)
	} ;
	if timings { println!("coeff = {}", coeff) ; }


        if timings { println!("BEFORE data: {} ms", now.elapsed().as_millis()); }

        let data =
            (&indptr[..dim as usize]).iter()
            .map(|ind| {
	        if debug { println!("indp[] = {}", ind) ; }
	        if (ind & z_indices).count_ones() % 2 == 1 {
		    -coeff
	        }
	        else {
		    coeff
	        }
            })
            .collect() ;

        if timings { println!("AFTER data: {} ms", now.elapsed().as_millis()); }

	Ok(UnsafeVectors {
            data,
            indices,
            indptr,
	})
    }
}

pub mod rowwise {
    use si_scale::helpers::{seconds, number_};
    use std::time::Instant;
    use num_complex::Complex64;
    use num_traits::Zero;
    use rayon::prelude::*;
    use sprs::{TriMatI};
    use std::cmp::min;
    use conv::prelude::* ;

    use super::UnsafeVectors ;

    pub fn make_params(members : &[crate::PauliSummand],
    ) -> Vec<(u64, u64, Complex64)> {
        members.iter()
            .map(|(p, coeff)| {
                let phase = p.phase() ;
                let coeff = *coeff ;
	        let coeff = match phase % 4 {
	            0 => Complex64::new(1.0, 0.0) * coeff,
	            1 => Complex64::new(0.0, -1.0) * coeff,
	            2 => Complex64::new(-1.0, 0.0) * coeff,
	            3 => Complex64::new(0.0, 1.0) * coeff,
	            _ => coeff // really should be assert!(false)
	        } ;
                (p.z_indices(), p.x_indices(), coeff)
            })
            .collect()
    }

    type RowContents = (Vec<u64>, Vec<Complex64>) ;

    pub fn make_row(params : &Vec<(u64, u64, Complex64)>, rowind : u64)
        -> RowContents
    {
        let mut v : Vec<(u64, Complex64)> = params.iter()
            .map(|(z_indices, x_indices, coeff)| {
                let colind = rowind ^ x_indices ;
                let coeff = *coeff ;
                let minus_coeff = -coeff ;
                let data = if (rowind & z_indices).count_ones() % 2 == 1 {
		    minus_coeff
	        }
	        else {
		    coeff
	        } ;
                (colind as u64, data)
            })
            .collect() ;
        v.sort_by(|a, b| a.0.cmp(&b.0)) ;
        let mut colw = Vec::new() ;
        let mut dataw = Vec::new() ;
        let (first, col, sum) = v.iter()
            .fold((true, 0, Complex64::zero()),
                  |(first, prevcol,runningsum),(colind,v)| {
                      if first {
                          (false, *colind, *v)
                      }
                      else if prevcol == *colind {
                          let runningsum = runningsum + *v ;
                          (false, prevcol, runningsum)
                      } else {
                          colw.push(prevcol) ;
                          dataw.push(runningsum) ;
                          (false, *colind, *v)
                      }
                  }) ;
        assert!(!first) ;
        colw.push(col as u64) ;
        dataw.push(sum) ;
        (colw, dataw)
    }

    pub fn make_trimat(members : &[crate::PauliSummand],
    ) -> TriMatI<Complex64,u64> {

        let num_qubits = members[0].0.num_qubits() ;
        let dim =  1 << num_qubits ;

        let params = make_params(members) ;

        let mut trimat : TriMatI<Complex64, u64> = TriMatI::new((dim, dim)) ;
        for rowind in 0..(dim as u64) {
            let rcont = make_row(&params, rowind) ;
            rcont.0.iter()
                .zip(rcont.1.iter())
                .for_each(|(colind, data)| {
                    trimat.add_triplet(rowind as usize, *colind as usize, *data)
                }) ;
        }
        trimat
    }

    pub fn make_unsafe_vectors(members : &[crate::PauliSummand],
    ) -> UnsafeVectors {

        let num_qubits = members[0].0.num_qubits() ;
        let dim =  1 << num_qubits ;

        let params = make_params(members) ;

        let accum_pairs : Vec<RowContents> =
            (0..(dim as u64)).into_par_iter()
             .map(|rowind| make_row(&params, rowind))
             .collect() ;

        let mut indptr = Vec::with_capacity(dim+1) ;
        let mut nnz : u64 = 0 ;
        for rowind in 0..(dim as u64) {
            let rc = &accum_pairs[rowind as usize] ;
            indptr.push(nnz) ;
            nnz += rc.0.len() as u64;
        }
        indptr.push(nnz) ;
        let mut indices = Vec::with_capacity(nnz as usize) ;
        let mut data = Vec::with_capacity(nnz as usize) ;
        accum_pairs.iter()
            .for_each(|(colv,datav)| {
                colv.iter().for_each(|colind| indices.push(*colind)) ;
		datav.iter().for_each(|v| data.push(*v)) ;
            }) ;
        UnsafeVectors {
            data,
            indices,
            indptr,
	}
    }

    pub fn make_unsafe_vectors_chunked(members : &[crate::PauliSummand],
                                       step : usize,
    ) -> UnsafeVectors {

        let debug = true ;
        let timings = true ;

        let now = Instant::now();

        if timings { println!("START make_unsafe_vectors_chunked: {}", seconds(now.elapsed().as_secs_f64())) ; }

        let num_qubits = members[0].0.num_qubits() ;
        let dim : usize =  1 << num_qubits ;

        let params = make_params(members) ;

        let chunks : Vec<(u64, u64)> =
            (0..(dim as u64))
            .into_iter()
            .step_by(step)
            .map(|n| (n,min(n + step as u64, dim as u64)))
            .collect();

        let chunked_vec : Vec<(Vec<u64>, Vec<RowContents>)> =
            chunks.par_iter()
            .map(|(lo,hi)| {
                let v_rc : Vec<RowContents> = (*lo..*hi).map(|rowind| make_row(&params, rowind)).collect() ;
		let v_nnz : Vec<u64> = v_rc.iter()
		    .map(|v| v.0.len() as u64).collect() ;
		let sum_nnz : u64 = v_nnz.iter().sum() ;
/*
		let mut indices = Vec::with_capacity(sum_nnz as usize) ;
		let mut data = Vec::with_capacity(sum_nnz as usize) ;
*/
		(v_nnz, v_rc)
            })
            .collect() ;

        if timings { println!("AFTER CHUNKS make_unsafe_vectors_chunked: {}", seconds(now.elapsed().as_secs_f64())) ; }

        let mut indptr = Vec::with_capacity(dim+1) ;
        let mut nnz : u64 = 0 ;
        for rowind in 0..(dim as u64) {
            let chunkind = rowind / step as u64 ;
            let chunkofs = rowind % step as u64 ;
	    let row_nnz = chunked_vec[chunkind as usize].0[chunkofs as usize] ;
            indptr.push(nnz) ;
            nnz += row_nnz as u64;
        }
        if debug { println!("EVENT nnz: {}", number_(f64::value_from(nnz).unwrap())) ; }
        indptr.push(nnz) ;
        if timings { println!("AFTER built indptr: {}", seconds(now.elapsed().as_secs_f64())) ; }
        let mut indices = Vec::with_capacity(nnz as usize) ;
        let mut data = Vec::with_capacity(nnz as usize) ;
        chunked_vec.iter()
            .for_each(|(_,vv)|
                      vv.iter()
                      .for_each(|(colv,datav)| {
			  colv.iter().for_each(|colind| indices.push(*colind)) ;
			  datav.iter().for_each(|v| data.push(*v)) ;
		      })
                      ) ;
        if debug { println!("EVENT indices.len()={} data.len()={} ",
			    number_(f64::value_from(indices.len()).unwrap()),
			    number_(f64::value_from(data.len()).unwrap())) ; }
        if timings { println!("END make_unsafe_vectors_chunked: {}", seconds(now.elapsed().as_secs_f64())) ; }
        UnsafeVectors {
            data,
            indices,
            indptr,
	}
    }
}
