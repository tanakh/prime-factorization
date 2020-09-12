use std::convert::TryInto;

use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use num_bigint::RandBigInt;
use prime_factorization::{factorization, gen_prime, is_prime};

fn millar_rabin_bench(c: &mut Criterion) {
    let mut rng = rand::thread_rng();

    let cases = [16, 32, 64, 128, 1024, 2048, 4096];

    for &bits in cases.iter() {
        c.bench_function(&format!("millar_rabin {}bit", bits), |b| {
            b.iter_batched(
                || rng.gen_biguint(bits),
                |n| is_prime(&n),
                BatchSize::LargeInput,
            )
        });
    }
}

fn gen_semiprime(bits: usize) -> u128 {
    let mut rng = rand::thread_rng();
    (gen_prime(bits, &mut rng) * gen_prime(bits, &mut rng))
        .try_into()
        .unwrap()
}

fn factorization_bench(c: &mut Criterion) {
    let cases = [8, 16, 24, 32, 40, 48];

    for &bits in cases.iter() {
        c.bench_function(&format!("factor semiprime {}bit", bits * 2), |b| {
            b.iter_batched(
                || gen_semiprime(bits),
                |n| factorization(n),
                BatchSize::LargeInput,
            )
        });
    }
}

criterion_group!(benches, millar_rabin_bench, factorization_bench);
criterion_main!(benches);
