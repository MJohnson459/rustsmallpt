use utility::*;

use nalgebra::core::Vector3;

pub type Vec3d = Vector3<f64>;

pub fn mult(vec1: &Vec3d, vec2: &Vec3d) -> Vec3d {
    Vec3d::new(
        vec1.x * vec2.x,
        vec1.y * vec2.y,
        vec1.z * vec2.z,
        )
}


pub fn from_rgb(hex: u32) -> Vec3d {
    Vec3d::new(
        ((hex & 0xFF0000) >> 16) as f64 / 255.0,
        ((hex & 0x00FF00) >> 8) as f64 / 255.0,
        ((hex & 0x0000FF) >> 0) as f64 / 255.0,
    )
}


pub fn vec_clamp(vec: &Vec3d) -> Vec3d {
    Vec3d::new(
        clamp(vec.x),
        clamp(vec.y),
        clamp(vec.z),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rbg() {
        let rgb = from_rgb(0xFF0000);
        let vec = Vec3d::new(1.0, 0.0, 0.0);
        assert_eq!(rgb, vec);
    }

    #[test]
    fn test_normalize() {
        let v: Vec3d = Vec3d {
            x: -0.0,
            y: 0.001017,
            z: -0.0,
        };
        let w: Vec3d = Vec3d {
            x: 0.0,
            y: -0.003064,
            z: 0.0,
        }.normalize();

        let u: Vec3d = v.normalize();

        assert!((u.length() - 1.0).abs() < 0.001);
        assert!((w.length() - 1.0).abs() < 0.001);
        println!("u: {:?}", u);
        println!("w: {:?}", w);

        let p: Vec3d = Vec3d {
            x: 0.0,
            y: -0.002908,
            z: 0.0,
        }.normalize();

        println!("p: {:?}", p);
        assert!((p.length() - 1.0).abs() < 0.001);

        let a = Vec3d {
            x: 0.003064,
            y: 0.0,
            z: -0.003064,
        }.normalize();
        let b = Vec3d {
            x: 0.003064,
            y: 0.0,
            z: -0.003064,
        }.normalize();
        println!("a: {:?}, a.length: {}", a, a.length());
        println!("b: {:?}, b.length: {}", b, b.length());
        assert!((a.length() - 1.0).abs() < 0.001);
        assert!((b.length() - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_cross() {
        let u: Vec3d = Vec3d {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        };
        let v: Vec3d = Vec3d {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        };

        assert_eq!(
            u.cross(v),
            Vec3d {
                x: 0.0,
                y: 0.0,
                z: -1.0,
            }
        );
        println!("u: {:?}", u);
        println!("v: {:?}", v);

        let w = Vec3d {
            x: 0.003064,
            y: 0.999991,
            z: 0.003064,
        };

        let u1 = Vec3d {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }.cross(w);
        let u2 = Vec3d {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }.cross(w);

        assert_eq!(
            u1,
            Vec3d {
                x: 0.003064,
                y: 0.0,
                z: -0.003064,
            }
        );
        assert_eq!(
            u2,
            Vec3d {
                x: 0.0,
                y: -0.003064,
                z: 0.999991,
            }
        );
    }

    /*#[test]
    fn bench_normalize() {
        let mut x = Vec3d::new(140.5,13.6,127.4);
        let time_start = time::precise_time_s();
        for _ in 0..10000000 {
            x.normalize();
        }
        let final_time = time::precise_time_s() - time_start;
        println!("final_time: {}", final_time);
        assert!(final_time < 0.1);
    }*/
}
