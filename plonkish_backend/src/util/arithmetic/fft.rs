//! Copy of https://github.com/privacy-scaling-explorations/halo2curves/blob/main/src/fft.rs.

use crate::util::{
    arithmetic::{Field, GroupOpsOwned, ScalarMulOwned},
    parallel::{join, num_threads},
    start_timer,
};

pub trait FftGroup<Scalar: Field>:
    Copy + Send + Sync + 'static + GroupOpsOwned + ScalarMulOwned<Scalar>
{
}

impl<T, Scalar> FftGroup<Scalar> for T
where
    Scalar: Field,
    T: Copy + Send + Sync + 'static + GroupOpsOwned + ScalarMulOwned<Scalar>,
{
}

pub fn radix2_fft<Scalar: Field, G: FftGroup<Scalar>>(a: &mut [G], omega: Scalar, log2_n: usize) {
    let _timer = start_timer(|| "fft");

    fn bitreverse(mut n: usize, l: usize) -> usize {
        let mut r = 0;
        for _ in 0..l {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    let log_num_threads = num_threads().ilog2() as usize;
    let n = a.len();
    assert_eq!(n, 1 << log2_n);

    for k in 0..n {
        let rk = bitreverse(k, log2_n);
        if k < rk {
            a.swap(rk, k);
        }
    }

    let twiddles: Vec<_> = (0..(n / 2))
        .scan(Scalar::ONE, |w, _| {
            let tw = *w;
            *w *= &omega;
            Some(tw)
        })
        .collect();

    if log2_n <= log_num_threads {
        let mut chunk = 2;
        let mut twiddle_chunk = n / 2;
        for _ in 0..log2_n {
            a.chunks_mut(chunk).for_each(|coeffs| {
                let (left, right) = coeffs.split_at_mut(chunk / 2);

                let (a, left) = left.split_at_mut(1);
                let (b, right) = right.split_at_mut(1);
                let t = b[0];
                b[0] = a[0];
                a[0] += &t;
                b[0] -= &t;

                left.iter_mut()
                    .zip(right.iter_mut())
                    .enumerate()
                    .for_each(|(i, (a, b))| {
                        let mut t = *b;
                        t *= &twiddles[(i + 1) * twiddle_chunk];
                        *b = *a;
                        *a += &t;
                        *b -= &t;
                    });
            });
            chunk *= 2;
            twiddle_chunk /= 2;
        }
    } else {
        recursive_butterfly_arithmetic(a, n, 1, &twiddles)
    }
}

fn recursive_butterfly_arithmetic<Scalar: Field, G: FftGroup<Scalar>>(
    a: &mut [G],
    n: usize,
    twiddle_chunk: usize,
    twiddles: &[Scalar],
) {
    if n == 2 {
        let t = a[1];
        a[1] = a[0];
        a[0] += &t;
        a[1] -= &t;
    } else {
        let (left, right) = a.split_at_mut(n / 2);
        join(
            || recursive_butterfly_arithmetic(left, n / 2, twiddle_chunk * 2, twiddles),
            || recursive_butterfly_arithmetic(right, n / 2, twiddle_chunk * 2, twiddles),
        );

        let (a, left) = left.split_at_mut(1);
        let (b, right) = right.split_at_mut(1);
        let t = b[0];
        b[0] = a[0];
        a[0] += &t;
        b[0] -= &t;

        left.iter_mut()
            .zip(right.iter_mut())
            .enumerate()
            .for_each(|(i, (a, b))| {
                let mut t = *b;
                t *= &twiddles[(i + 1) * twiddle_chunk];
                *b = *a;
                *a += &t;
                *b -= &t;
            });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::arithmetic::{root_of_unity, root_of_unity_inv};
    use halo2_curves::bn256::Fr;
    use rand::thread_rng;

    #[test]
    fn test_radix2_fft() {
        let n: u64 = 8;
        let n_inv = Fr::from(n).invert().unwrap();
        assert!(n_inv * Fr::from(n) == Fr::ONE);
        let log_n = 3; // 2^3 = 8
                       // generate test data
        let mut rng = thread_rng();
        let mut data: Vec<Fr> = (0..n).map(|_| Fr::random(&mut rng)).collect();
        let original_data = data.clone();

        // calculate appropriate unit root
        let omega = root_of_unity(n as usize);
        let omega_inv = root_of_unity_inv(n as usize);

        // execute FFT
        radix2_fft(&mut data, omega, log_n);
        println!("fft result: {:?}", data);

        // execute IFFT to convert fft result back to original data
        radix2_fft(&mut data, omega_inv, log_n);
        // divide by n to get ifft result
        data.iter_mut().for_each(|data| *data *= n_inv);
        println!("ifft result: {:?}", data);

        assert_eq!(data, original_data);

        println!("radix2_fft test passed");
    }
}
