use halo2_curves::bn256::Fr;
use std::ops::{Add, Mul};

type Scalar = Fr;

pub struct SOSTable {
    pub(crate) l: usize,
    pub(crate) c: usize,
    pub(crate) k: usize,
    pub(crate) alpha: usize,
    pub(crate) tables: Vec<Vec<i64>>,
}

pub struct Params {
    pub(crate) table: SOSTable,
}

impl SOSTable {
    pub fn new(l: usize, c: usize, k: usize, tables: Vec<Vec<i64>>) -> Self {
        assert_eq!(tables.len(), k * c);
        Self {
            l,
            c,
            k,
            alpha: k * c,
            tables,
        }
    }

    // T could be a Scalar
    pub fn g_func<T>(&self, r: &[T]) -> T
    where
        T: Clone + Add<Output = T> + Mul<Output = T> + From<i64>,
    {
        /*
        The g function.
            - r could be a list of Scalar or poly.Polynomial or mle_poly.polynomial.
            - ret would be the same type as r_1.
        Here we have g(r_1, ..., r_c) = 1 + r_1 + r_2 * 2**l + ... + r_c * 2**(l*(c-1))
        that represent the table [1, 2, ..., 2**(l*c)].
        [Need to override this function.]
         */
        assert_eq!(r.len(), self.alpha);

        let mut ret = T::from(1);
        let mut mul = T::from(1);
        let base = T::from(2i64.pow(self.l as u32));

        for ri in r {
            ret = ri.clone() * mul.clone() + ret;
            mul = mul * base.clone();
        }
        return ret;
    }

    pub fn get_index(&self, value: i64) -> Vec<i64> {
        /*
        Return the index of each subtable for any element in the table.
        With the SOS structure, we should be able to get the indexes
        without iterating through the whole table.
        Here we return the indexes for the element in the table [1, 2, ..., 2**(l*c)]
        [Need to override this function]
        */
        let mut val = value - 1;
        let mut index = Vec::new();
        for _ in 0..self.c {
            index.push(val % (1 << self.l));
            val /= 1 << self.l;
        }
        return index;
    }
}

impl Params {
    pub fn new(table: SOSTable) -> Self {
        Self { table }
    }
}
