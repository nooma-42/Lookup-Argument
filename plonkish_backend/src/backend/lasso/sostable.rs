use crate::poly::multilinear::MultilinearPolynomial;
use halo2_curves::bn256::Fr;

type Scalar = Fr;

pub struct SOSTable {
    pub(crate) l: usize,
    pub(crate) c: usize,
    pub(crate) k: usize,
    pub(crate) alpha: usize,
    pub(crate) tables: Vec<Vec<u64>>,
}

pub struct Params {
    pub(crate) table: SOSTable,
}

impl SOSTable {
    pub fn new(l: usize, c: usize, k: usize, tables: Vec<Vec<u64>>) -> Self {
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
    pub fn g_func(&self, r: &[MultilinearPolynomial<Scalar>]) -> MultilinearPolynomial<Scalar> {
        /*
        The g function.
        Here we have g(r_1, ..., r_c) = 1 + r_1 + r_2 * 2**l + ... + r_c * 2**(l*(c-1))
        that represent the table [1, 2, ..., 2**(l*c)].
        [Need to override this function.]
         */

        // doing fast pow for polynomial
        assert_eq!(r.len(), self.alpha);

        let mut ret = MultilinearPolynomial::new(vec![Scalar::one()]);
        let mut pow = Scalar::from(1 << self.l); // 2^l as Scalar
        let base = Scalar::from(1 << self.l); // 2^l as Scalar
        for poly in r.iter() {
            ret = &(poly * pow) + &ret;
            pow = &pow * &base;
        }
        ret
    }

    pub fn get_index(&self, value: u64) -> Vec<u64> {
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
