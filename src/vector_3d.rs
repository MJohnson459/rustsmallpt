use std::f64;

use std::ops::{Add,Sub,Mul,Rem};


// ------- VEC ------------
#[derive(Copy, Clone, Debug)]
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

	pub fn normalise(&mut self) -> Vec3d {
		let nx = self.x.powi(2);
		let ny = self.y.powi(2);
		let nz = self.z.powi(2);
		*self = *self * (1.0/(nx + ny + nz).sqrt());
		*self
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
		(self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
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

impl Rem for Vec3d {
	type Output = Vec3d;
	fn rem(self, other: Vec3d) -> Vec3d {
		Vec3d {
			x: self.y * other.z - self.z*other.y,
			y: self.z * other.x - self.x*other.z,
			z: self.x * other.y - self.y*other.x,
		}
	}
}

impl<'a,'b> Rem<&'b Vec3d> for &'a Vec3d {
	type Output = Vec3d;
	fn rem(self, other: &'b Vec3d) -> Vec3d {
		Vec3d {
			x: self.y * other.z - self.z*other.y,
			y: self.z * other.x - self.x*other.z,
			z: self.x * other.y - self.y*other.x,
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
