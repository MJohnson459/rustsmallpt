use std::f64;

use std::ops::{Add,Sub,Mul,AddAssign};


// ------- VEC ------------
#[derive(Default, Copy, Clone, Debug)]
pub struct Vec3d {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3d {
    pub fn new(x: f64, y: f64, z: f64) -> Vec3d {
        Vec3d {x: x, y: y, z: z}
    }

    pub fn zeros() ->Vec3d {
        Vec3d {x: 0.0, y: 0.0, z: 0.0}
    }

    pub fn normalise(&self) -> Vec3d {
        let nx = self.x * self.x;
        let ny = self.y * self.y;
        let nz = self.z * self.z;
        *self * (1.0/(nx + ny + nz).sqrt())
    }

    pub fn dot(&self, other: &Vec3d) -> f64 {
        self.x*other.x + self.y*other.y + self.z*other.z
    }

    pub fn cross(self, other: Vec3d) -> Vec3d {
        Vec3d {
            x: self.y*other.z - self.z*other.y,
            y: self.z*other.x - self.x*other.z,
            z: self.x*other.y - self.y*other.x,
        }
    }

    pub fn length(self) -> f64 {
        (self.x*self.x + self.y*self.y + self.z*self.z).sqrt()
    }

    pub fn mult(self, other: Vec3d) -> Vec3d {
        Vec3d {
            x: self.x * other.x,
            y: self.y * other.y,
            z: self.z * other.z,
        }

    }
}

impl Add for Vec3d {
    type Output = Vec3d;
    fn add(self, other: Vec3d) -> Vec3d {
        Vec3d {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl AddAssign for Vec3d {
    fn add_assign(&mut self, other: Vec3d) {
        *self = Vec3d {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        };
    }
}

impl<'a,'b> Sub<&'b Vec3d> for &'a Vec3d {
    type Output = Vec3d;
    fn sub(self, other: &'b Vec3d) -> Vec3d {
        Vec3d {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Sub for Vec3d {
    type Output = Vec3d;
    fn sub(self, other: Vec3d) -> Vec3d {
        Vec3d {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a,'b> Mul<&'b Vec3d> for &'a Vec3d {
    type Output = Vec3d;
    fn mul(self, other: &'b Vec3d) -> Vec3d {
        Vec3d {
            x: self.x * other.x,
            y: self.y * other.y,
            z: self.z * other.z,
        }
    }
}

impl Mul<Vec3d> for Vec3d {
    type Output = Vec3d;
    fn mul(self, other: Vec3d) -> Vec3d {
        Vec3d {
            x: self.x * other.x,
            y: self.y * other.y,
            z: self.z * other.z,
        }
    }
}

impl Mul<f64> for Vec3d {
    type Output = Vec3d;
    fn mul(self, other: f64) -> Vec3d {
        Vec3d {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl PartialEq for Vec3d {
    fn eq(&self, other: &Vec3d) -> bool {
        ((self.x - other.x).abs() < 0.0001) &&
        ((self.y - other.y).abs() < 0.0001) &&
        ((self.z - other.z).abs() < 0.0001)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalise() {
        let v: Vec3d = Vec3d{x: -0.0, y: 0.001017, z: -0.0};
        let w: Vec3d = Vec3d{x: 0.0, y: -0.003064, z: 0.0}.normalise();

        let u: Vec3d = v.normalise();

        assert!((u.length() - 1.0).abs() < 0.001);
        assert!((w.length() - 1.0).abs() < 0.001);
        println!("u: {:?}", u);
        println!("w: {:?}", w);

        let p: Vec3d = Vec3d{x: 0.0, y: -0.002908, z: 0.0}.normalise();

        println!("p: {:?}", p);
        assert!((p.length() - 1.0).abs() < 0.001);

        let a = Vec3d{x: 0.003064, y: 0.0, z: -0.003064}.normalise();
        let b = Vec3d{x: 0.003064, y: 0.0, z: -0.003064}.normalise();
        println!("a: {:?}, a.length: {}", a, a.length());
        println!("b: {:?}, b.length: {}", b, b.length());
        assert!((a.length() - 1.0).abs() < 0.001);
        assert!((b.length() - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_cross() {
        let u: Vec3d = Vec3d{x: 0.0, y: 1.0, z: 0.0};
        let v: Vec3d = Vec3d{x: 1.0, y: 0.0, z: 0.0};

        assert_eq!(u.cross(v), Vec3d{x: 0.0, y: 0.0, z: -1.0});
        println!("u: {:?}", u);
        println!("v: {:?}", v);

        let w = Vec3d { x: 0.003064, y: 0.999991, z: 0.003064 };

        let u1 = Vec3d{x:0.0,y:1.0,z:0.0}.cross(w);
        let u2 = Vec3d{x:1.0,y:0.0,z:0.0}.cross(w);

        assert_eq!(u1, Vec3d{x: 0.003064, y: 0.0, z: -0.003064});
        assert_eq!(u2, Vec3d{x: 0.0, y: -0.003064, z: 0.999991});
    }


    /*#[test]
    fn bench_normalise() {
        let mut x = Vec3d::new(140.5,13.6,127.4);
        let time_start = time::precise_time_s();
        for _ in 0..10000000 {
            x.normalise();
        }
        let final_time = time::precise_time_s() - time_start;
        println!("final_time: {}", final_time);
        assert!(final_time < 0.1);
    }*/
}
