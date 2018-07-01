use utility::*;

use cgmath::Vector3;

pub type Vec3d = Vector3<f64>;


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
