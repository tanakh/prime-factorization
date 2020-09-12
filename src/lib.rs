use num_bigint::{BigUint, RandBigInt};
use rand::Rng;

uint::construct_uint! {
    struct U256(4);
}

pub fn factorization(n: u128) -> Vec<u128> {
    if is_prime(&n.into()) {
        return vec![n];
    }

    let mut ret = vec![];
    let mut n = n;

    while !is_prime(&n.into()) {
        let factor = pollard_rho(n);
        assert!(n > factor);
        assert!(n % factor == 0);
        ret.append(&mut factorization(factor));
        n /= factor;
    }

    if n > 1 {
        ret.push(n);
    }

    ret.sort();
    ret
}

pub fn is_prime(n: &BigUint) -> bool {
    if n == &0u32.into() || n == &1u32.into() {
        false
    } else if n == &2u32.into() || n == &3u32.into() {
        true
    } else {
        miller_rabin_test(n, 100) == MillerRabinResult::ProbablyPrime
    }
}

#[derive(PartialEq, Eq, Debug)]
enum MillerRabinResult {
    Composite,
    ProbablyPrime,
}

fn miller_rabin_test(n: &BigUint, k: usize) -> MillerRabinResult {
    let t: BigUint = n - 1_u32;
    let r = t.trailing_zeros().unwrap();
    let d = t.clone() >> r;
    let one = BigUint::from(1u32);
    let two = BigUint::from(2u32);

    let mut rng = rand::thread_rng();
    'outer: for _ in 0..k {
        let a: BigUint = rng.gen_biguint_range(&two, &t);
        let mut x = a.modpow(&d, &n);

        if x == one || x == t {
            continue;
        }
        for _ in 1..r {
            x = x.modpow(&two, &n);
            if x == t {
                continue 'outer;
            }
        }
        return MillerRabinResult::Composite;
    }
    MillerRabinResult::ProbablyPrime
}

pub fn gen_prime(bits: usize, rng: &mut impl Rng) -> BigUint {
    loop {
        let n = rng.gen_biguint(bits as _);
        if is_prime(&n) {
            break n;
        }
    }
}

pub fn pollard_rho(n: u128) -> u128 {
    loop {
        if let Some(ret) = pollard_rho_once(n, rand::thread_rng().gen_range(2, n)) {
            return ret;
        }
    }
}

fn pollard_rho_once(n: u128, mut x: u128) -> Option<u128> {
    for cycle in 1.. {
        let y = x;
        for _ in 0..1 << cycle {
            x = mul_mod(x, x, n);
            x = (x + 1) % n;
            let d = if x >= y { x - y } else { y - x };
            let factor = gcd(d, n);
            if factor > 1 {
                return if factor != n { Some(factor) } else { None };
            }
        }
    }
    unreachable!()
}

fn gcd(a: u128, b: u128) -> u128 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn mul_mod(a: u128, b: u128, m: u128) -> u128 {
    (U256::from(a) * U256::from(b) % U256::from(m)).as_u128()
}

#[cfg(test)]
mod tests {
    use crate::*;
    use rand::Rng;
    use MillerRabinResult::*;

    #[test]
    fn miller_rabin() {
        assert_eq!(miller_rabin_test(&4u32.into(), 100), Composite);
        assert_eq!(miller_rabin_test(&5u32.into(), 100), ProbablyPrime);
        assert_eq!(miller_rabin_test(&6u32.into(), 100), Composite);
        assert_eq!(miller_rabin_test(&7u32.into(), 100), ProbablyPrime);
        assert_eq!(miller_rabin_test(&8u32.into(), 100), Composite);
        assert_eq!(miller_rabin_test(&9u32.into(), 100), Composite);
        assert_eq!(miller_rabin_test(&10u32.into(), 100), Composite);

        let primes = [
            241393502644931236824083437316691947053_u128,
            203851774287909279562700874830405601907,
            288998372467163468683539125188698897449,
            189822170223895762265064303540250435637,
            328089814025503192461392250603431433007,
            296924869529206030921754222802047897761,
            239796305212141579879826223504530404303,
            183793893507075905518837979761809022117,
        ];

        for &p in primes.iter() {
            assert_eq!(miller_rabin_test(&p.into(), 100), ProbablyPrime);
        }

        let composites = [
            123818405246410390178707765938870261443_u128,
            202383583825519051662160491314347582777,
            167210971201646539721622042892407080737,
            117531937261998017622502744988913176479,
            212783729721844122818518594528126775039,
            305792237550134519087322660700463788397,
            93694064147469858944238195598447143569,
            132633955283692567229336700018758707763,
        ];

        for &c in composites.iter() {
            assert_eq!(miller_rabin_test(&c.into(), 100), Composite);
        }
    }

    #[test]
    fn factorization_test() {
        assert_eq!(factorization(2), [2]);
        assert_eq!(factorization(3), [3]);
        assert_eq!(factorization(4), [2, 2]);
        assert_eq!(factorization(5), [5]);
        assert_eq!(factorization(6), [2, 3]);
        assert_eq!(factorization(7), [7]);
        assert_eq!(factorization(8), [2, 2, 2]);
        assert_eq!(factorization(9), [3, 3]);
        assert_eq!(factorization(10), [2, 5]);

        let test = |x| {
            let fs = factorization(x);
            assert!(fs.iter().all(|&f| is_prime(&f.into())));
            assert_eq!(fs.iter().product::<u128>(), x);
        };

        for x in 2..1000 {
            test(x);
        }

        for _ in 0..100 {
            let x = rand::thread_rng().gen_range(2, u64::MAX as u128);
            test(x);
        }
    }
}
