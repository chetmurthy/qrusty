use num_complex::Complex64;
use std::time::Instant;

pub fn rust_make_data(z: Vec<bool>,
                      x: Vec<bool>,
                      coeff: Complex64,
		      phase: i64,
		      group_phase: bool
                ) -> std::result::Result<(Vec<Complex64>, Vec<u64>, Vec<u64>), &'static str> {
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
	    let mut dotprod = 0 ;
	    for i in 0..num_qubits {
		if x[i] && z[i] {
		    dotprod += 1
		}
	    }
	    if debug { println!("2: dotprod={}", dotprod) ; }
	    mut_phase += dotprod ;
	    mut_phase = mut_phase % 4 ;
	}

	if debug { println!("2: mut_phase={}", mut_phase) ; }

	let dim =  1 << num_qubits ;
	let mut twos_array = Vec::<u64>::new() ;
	for i in 0..num_qubits {
	    twos_array.push(1 << i) ;
	}
	if debug { println!("3: twos_array={:?}", twos_array) ; }

	let mut x_indices = 0 ;
	for i in 0..num_qubits {
	    if x[i] {
		x_indices += twos_array[i] ;
	    }
	}
	let mut z_indices = 0 ;
	for i in 0..num_qubits {
	    if z[i] {
		z_indices += twos_array[i] ;
	    }
	}

	if debug { println!("4: x_indices={} z_indices={}", x_indices, z_indices) ; }


        if timings { println!("BEFORE indptr: {} ms", now.elapsed().as_millis()); }


	let mut indptr = vec![0 as u64;dim+1] ;
	for i in 0..(dim+1) {
	    indptr[i] = i as u64;
	}

        if timings { println!("BEFORE indices: {} ms", now.elapsed().as_millis()); }
	let mut indices = vec![0 as u64;indptr.len()] ;

	for i in 0..indptr.len() {
	    indices[i] = indptr[i] ^ x_indices ;
	}
	let coeff = match phase % 4 {
	    0 => Complex64::new(1.0, 0.0) * coeff,
	    1 => Complex64::new(0.0, -1.0) * coeff,
	    2 => Complex64::new(-1.0, 0.0) * coeff,
	    3 => Complex64::new(0.0, 1.0) * coeff,
	    _ => coeff // really should be assert!(false)
	} ;
	if timings { println!("coeff = {}", coeff) ; }


        if timings { println!("BEFORE data: {} ms", now.elapsed().as_millis()); }

	let mut data = Vec::new() ;
	for indp in indptr.iter() {
	    if debug { println!("indp[] = {}", indp) ; }
	    if (indp & z_indices).count_ones() % 2 == 1 {
		data.push(-coeff) ;
	    }
	    else {
		data.push(coeff) ;
	    }
	}
        if timings { println!("AFTER data: {} ms", now.elapsed().as_millis()); }

	Ok((
            data,
            indices,
            indptr,
	))
    }
}
