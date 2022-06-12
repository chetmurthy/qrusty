use std::boxed::Box ;

pub struct BinaryTreeFold<T> {
    stk : Vec<(usize, T)>,
    op : fn(l : T, r: T) -> T
}

/* a BinaryTreeFold is a stack that keeps track of where we are in the
 * binary-tree-based folding process.
 *
 * The operations are:
 *
 * begin: given a value of type T and a binary operator over T, create
 * the BTF
 *
 * add: add a new T to the BTF, invoking the binary operator, perhaps
 * more than once
 *
 * end: extract the final value, again perhaps invoking the binary
 * operator more than once
 */

impl<T> BinaryTreeFold<T> {
    pub fn begin(init : T, op : fn(l: T,r: T) -> T) -> BinaryTreeFold<T> {
        let stk = vec![ (1, init) ] ;
        BinaryTreeFold { stk, op }
    }
/*
            def addp(self, p):
                if not self.stack:
                    self.stack.append(p)
                    return
                q = self.stack[-1]
                if q[0] == p[0]:
                    height = 1 + q[0]
                    newv = self.f(height, q[1], p[1])
                    self.stack.pop()
                    self.addp((height, newv))
                else:
                    self.stack.append(p)

            def add(self, v):
                self.addp((1,v))

            def flush(self):
                if not self.stack:
                    raise Exception("flush: stack was empty")
                v = self.stack.pop()[1]
                while self.stack:
                    q = self.stack.pop()
                    v = self.f(q[0], q[1], v)
                return v

*/
    /* invariant: the stack is always nonempty */
    fn addp(&mut self, v : (usize, T)) {
        if self.stk.len() == 0 {
            self.stk.push(v) ;
        }
        else {
            let q = self.stk.last().unwrap() ;
            if q.0 == v.0 {
                let q = self.stk.pop().unwrap() ;
                let newv = (self.op)(q.1, v.1) ;
                self.addp((q.0 + 1, newv)) ;
            } else {
                self.stk.push(v) ;
            }
        }
    }

    pub fn add(&mut self, v : T) {
        self.addp((1,v))
    }

    pub fn end(&mut self) -> T {
        let mut v = self.stk.pop().unwrap().1 ;
        while self.stk.len() > 0 {
            let w = self.stk.pop().unwrap() ;
            v = (self.op)(v, w.1) ;
        }
        drop(self) ;
        v
    }
}

pub mod list {
    use std::fmt;

    #[derive(Debug)]
    pub enum ListNode<T> {
        Cons(List<T>, List<T>),
        Atom(T)
    }
    type List<T> = Box<ListNode<T>> ;


    pub fn cons<T>(l : List<T>, r : List<T>) -> List<T> {
        Box::new(ListNode::Cons(l,r))
    }

    pub fn atom<T>(a : T) -> List<T> {
        Box::new(ListNode::Atom(a))
    }

    impl<T: fmt::Display> fmt::Display for ListNode<T> {
        // This trait requires `fmt` with this exact signature.
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            match self {
                ListNode::Atom(t) => write!(f, "{}", t),
                ListNode::Cons(l,r) =>
                    write!(f, "({} {})", l,  r)
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use crate::util::BinaryTreeFold ;
    use crate::util::list::* ;

    #[test]
    fn plus() {
        let mut bt = BinaryTreeFold::begin(1, |x,y| x+y) ;
        (2..10).for_each(|n| bt.add(n)) ;

        assert_eq!((1..10).fold(0,|a,b| a+b), bt.end()) ;
    }

    #[test]
    fn test_cons() {
        let mut bt = BinaryTreeFold::begin(atom(0), cons) ;
        (1..8).for_each(|n| bt.add(atom(n))) ;
        let rv = bt.end() ;
        assert_eq!("(((0 1) (2 3)) ((4 5) (6 7)))",
                   (*rv).to_string()) ;
    }
}
