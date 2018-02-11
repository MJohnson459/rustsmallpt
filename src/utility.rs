extern crate time;

// pub fn float_eq(a: f64, b:f64) -> bool {
//     (a - b).abs() < 0.0001
// }

pub fn clamp(x: f64) -> f64 {
    if x < 0.0 {
        0.0
    } else if x > 1.0 {
        1.0
    } else {
        x
    }
}

// pub fn to_int(x: f64) -> i32 {
//     (clamp(x).powf(1.0/2.2)*255.0+0.5) as i32
// }

pub fn to_u8(x: f64) -> u8 {
    (gamma_correct(clamp(x)) * 255.0 + 0.5) as u8
}

fn gamma_correct(x: f64) -> f64 {
    x.powf(1.0 / 2.2)
}

pub fn format_time(seconds: f64) -> String {
    let hours: f64 = (seconds / 3600.0).floor();
    let mins: f64 = ((seconds - hours * 3600.0) / 60.0).floor();
    let secs: f64 = seconds - mins * 60.0;
    format!("{:02.0}:{:02.0}:{:02.0}", hours, mins, secs)
}

#[cfg(all(feature = "unstable", test))]
mod bench {
    extern crate test;

    use super::*;
    use self::test::Bencher;

    #[bench]
    fn bench_to_u8(b: &mut Bencher) {
        b.iter(|| {
            let n = test::black_box(0.65473);
            to_u8(n)
        });
    }

    #[bench]
    fn bench_gamma_correct(b: &mut Bencher) {
        b.iter(|| {
            let n = test::black_box(0.65473);
            gamma_correct(n)
        });
    }

}
