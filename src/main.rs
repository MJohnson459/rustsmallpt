extern crate rand;

use std::f64;

use std::ops::{Add,Sub,Mul,Rem};
use std::io::BufWriter;
use std::fs::{OpenOptions};
use std::io::Write;


// ------- VEC ------------
#[derive(Copy, Clone, Debug)]
pub struct Vec3d {
	x: f64,
	y: f64,
	z: f64,
}

impl Vec3d {
	fn normalise(self) -> Vec3d {
		let nx = self.x.powi(2);
		let ny = self.y.powi(2);
		let nz = self.z.powi(2);
		self * (1.0/(nx + ny + nz).sqrt())
	}

	fn dot(&self, other: &Vec3d) -> f64 {
		self.x*other.x + self.y*other.y + self.z*other.z
	}

	fn cross(self, other: Vec3d) -> Vec3d {
		Vec3d {
			x: self.y*other.z - self.z*other.y,
			y: self.z*other.x - self.x*other.z,
			z: self.x*other.y - self.y*other.x,
		}
	}

	fn mult(self, other: Vec3d) -> Vec3d {
		Vec3d {
			x: self.x * other.x,
			y: self.y * other.y,
			z: self.x * other.x,
		}

	}
}

impl Add for Vec3d {
	type Output = Vec3d;
	fn add(self, other: Vec3d) -> Vec3d {
		Vec3d {
			x: self.x + other.x,
			y: self.y + other.y,
			z: self.x + other.x,
		}
	}
}

impl<'a,'b> Sub<&'b Vec3d> for &'a Vec3d {
	type Output = Vec3d;
	fn sub(self, other: &'b Vec3d) -> Vec3d {
		Vec3d {
			x: self.x - other.x,
			y: self.y - other.y,
			z: self.x - other.x,
		}
	}
}

impl Sub for Vec3d {
	type Output = Vec3d;
	fn sub(self, other: Vec3d) -> Vec3d {
		Vec3d {
			x: self.x - other.x,
			y: self.y - other.y,
			z: self.x - other.x,
		}
	}
}

impl<'a,'b> Mul<&'b Vec3d> for &'a Vec3d {
	type Output = Vec3d;
	fn mul(self, other: &'b Vec3d) -> Vec3d {
		Vec3d {
			x: self.x * other.x,
			y: self.y * other.y,
			z: self.x * other.x,
		}
	}
}

impl Mul<f64> for Vec3d {
	type Output = Vec3d;
	fn mul(self, other: f64) -> Vec3d {
		Vec3d {
			x: self.x * other,
			y: self.y * other,
			z: self.x * other,
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

// ----------- Ray --------------

pub struct Ray {
	o: Vec3d,
	d: Vec3d,
}

// ----------
#[derive(Copy, Clone)]
enum ReflectType {
	DIFF,
	SPEC,
	REFR
}

// ------
#[derive(Copy, Clone)]
pub struct Sphere {
	radius: f64,
	position: Vec3d,
	emission: Vec3d,
	color: Vec3d,
	reflection: ReflectType
}

impl Sphere {
	fn intersect(&self, ray: &Ray) -> f64 {
		let op: Vec3d = &self.position - &ray.o;
		let eps: f64 = 1e-4;
		let b: f64 = op.dot(&ray.d);
		let mut det: f64 = b*b - op.dot(&op) + self.radius*self.radius;
		if det < 0.0 {
			//println!("det < 0.0 : {}", det);
			return 0.0
		} else {
			//println!("det >= 0.0 : {}", det);
			det = det.sqrt();
		}
		let t1 = b - det;
		let t2 = b + det;
		if t1 > eps {
			return t1;
		} else if t2 > eps {
			return t2;
		} else {
			//println!("t1 <= eps && t2 <= eps - t1: {}, t2: {}, eps: {}",t1,t2,eps );
			return 0.0;
		}
	}
}

// ---------------------

pub fn clamp(x: f64) -> f64 {
	if x < 0.0 {
		return 0.0;
	} else 	if x > 1.0 {
		return 0.0;
	}
	x
}

pub fn to_int(x: f64) -> i32 {
	(clamp(x).powf(1.0/2.2)*255.0+0.5) as i32
}

pub fn intersect(spheres: &mut Vec<Sphere>, ray: &Ray) -> (bool, f64, i32) {
	let n = spheres.len();
	let mut d: f64;
	let mut t: f64 = 1e20;
	let mut id: i32 = 0;
	let inf = t.clone();

	for i in (0..n).rev() {
		d = spheres[i].intersect(ray);
		if d != 0.0 && d < t {
			//println!("something intersected sphere i: {}", i);
			t = d;
			id = i as i32;
		}
	}
	//println!("intersects - t: {}, inf: {}", t, inf);
	(t < inf, t, id)
}

pub fn radiance(spheres: &mut Vec<Sphere>, ray: &Ray, mut depth: i32) -> Vec3d {
	println!("radiance - depth: {},  ray dir  x: {}, y: {}, z: {}", depth, ray.d.x, ray.d.y, ray.d.z);
	let (intersects, t , id) = intersect(spheres, ray);
	if !intersects {  // if miss, return black
		//println!("No Intersect!");
		return Vec3d{x:0.0,y:0.0,z:0.0};
	}

	let obj: Sphere = spheres[id as usize].clone();  // the hit object

	let x: Vec3d = ray.o + ray.d*t;
	let n: Vec3d = (&x - &obj.position).normalise();
	let nl: Vec3d;
	if n.dot(&ray.d) < 0.0 {
		nl = n;
	} else {
		nl = n*-1.0;
	}
	let mut f: Vec3d = obj.color.clone();

	let p: f64;
	if f.x > f.y && f.x>f.z {
		p = f.x;
	} else if f.y > f.z {
		p = f.y;
	} else {
		p = f.z;
	}

	depth = depth + 1;
	if depth > 5 {
		if rand::random::<f64>() < p {
			f = f*(1.0/p);
		} else {
			// This might be at the wrong if statement
			return obj.emission;
		}
	}

	match obj.reflection {
		ReflectType::DIFF => {
			let r1: f64 = 2.0*f64::consts::PI*0.5; //rand::random::<f64>();
			let r2: f64 = 0.5; //rand::random();
			let r2s = r2.sqrt();

			let w: Vec3d = nl.clone();
			let u: Vec3d;
			if w.x.abs() > 0.1 {
				u = (&Vec3d{x:0.0,y:1.0,z:0.0}%&w).normalise();
			} else {
				u = (&Vec3d{x:1.0,y:0.0,z:0.0}%&w).normalise();
			}

			let v: Vec3d = &w % &u;

			let d: Vec3d = (u*r1.cos()*r2s + v*r1.sin()*r2s + w*(1.0-r2).sqrt()).normalise();
			println!("obj.emission DIFF u: {:?}, v: {:?}, w: {:?}, r1: {:?}, r2s: {:?}", u, v, w, r1, r2s);
			let fmu = f.mult(radiance(spheres, &Ray{o:x,d:d}, depth));
			//println!("obj.emission DIFF - fmu.x: {}, fmu.y: {}, fmu.z: {}, depth: {}", fmu.x, fmu.y, fmu.z, depth);
			//println!("obj.emission.x: {}, obj.emission.y: {}, obj.emission.z: {}", obj.emission.x, obj.emission.y, obj.emission.z);
			return obj.emission + fmu;
		},
		ReflectType::SPEC => {
			//println!("obj.emission SPEC, depth: {}", depth);
			let d: Vec3d = ray.d-(n*2.0*n.dot(&ray.d));
			return obj.emission + f.mult(radiance(spheres, &Ray{o:x,d:d}, depth));
		},
		_ => {}

	}

	let refl_ray: Ray = Ray{o:x, d:ray.d-n*2.0*n.dot(&ray.d)}; // Ideal dielectric REFRACTION
	let into: bool = n.dot(&nl) > 0.0;

	let nc: f64 = 1.0;
	let nt: f64 = 1.5;
	let nnt: f64;
	if into { // Ray from outside going in?
		nnt = nc/nt;
	} else {
		nnt = nt/nc;
	}
	let ddn: f64 =ray.d.dot(&nl);
	let cos2t: f64 = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
	if cos2t < 0.0 { // Total internal reflection
		//println!("obj.emission cos2t < 0.0");
		return obj.emission + f.mult(radiance(spheres, &refl_ray,depth));
	}

	let tdir: Vec3d;

	if into {
		tdir = (ray.d*nnt - n*(1.0*(ddn*nnt+cos2t.sqrt()))).normalise();
	} else {
		tdir = (ray.d*nnt - n*(-1.0*(ddn*nnt+cos2t.sqrt()))).normalise();
	}

	let a: f64 = nt- nc;
	let b: f64 = nt + nc;
	let r0: f64 = a*a/(b*b);
	let c: f64;

	if into {
		c = 1.0 + ddn;
	} else {
		c = 1.0 - tdir.dot(&n);
	}

	let re: f64 = r0+(1.0-r0)*c.powi(5);
	let tr: f64 = 1.0 - re;
	let p: f64 = 0.25+0.5*re;
	let rp: f64 = re/p;
	let tp: f64 = tr/(1.0-p);

	if depth > 2 {
		if rand::random::<f64>() < p {
			//println!("rand::random::<f64>() < p - depth: {}", depth);
			return obj.emission + f.mult(radiance(spheres, &refl_ray,depth)*rp);
		} else {
			//println!("rand::random::<f64>() >= p - depth: {}", depth);
			return obj.emission + f.mult(radiance(spheres, &Ray{o:x,d:tdir},depth)*tp);
		}
	} else {
		//println!("reflect ray: {}", depth);
		return obj.emission + f.mult(radiance(spheres, &refl_ray,depth)*re+radiance(spheres, &Ray{o:x,d:tdir},depth)*tr);
	}
}

fn main() {
	println!("Hello, world!");

	let mut spheres = Vec::new();
	spheres.push(Sphere {radius:1e5, position: Vec3d{x:1e5+1.0,y:40.8,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.25,z:0.25}, reflection: ReflectType::DIFF}); // left
	spheres.push(Sphere {radius:1e5, position: Vec3d{x:-1e5+99.0,y:40.8,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.25,y:0.25,z:0.75}, reflection: ReflectType::DIFF}); // right
	spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:40.8,z:1e5}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // back
	spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:40.8,z:-1e5+170.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // front
	spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:1e5,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // bottom
	spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:-1e5+81.6,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // Top
	spheres.push(Sphere {radius:16.5, position: Vec3d{x:27.0,y:16.5,z:47.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.999, reflection: ReflectType::SPEC}); // mirror
	spheres.push(Sphere {radius:16.5, position: Vec3d{x:73.0,y:16.5,z:78.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.999, reflection: ReflectType::REFR}); // glass
	spheres.push(Sphere {radius:600.0, position: Vec3d{x:50.0,y:681.6-0.27,z:81.6}, emission: Vec3d{x:12.0,y:12.0,z:12.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // light

	let width: i32 = 100;
	let height: i32 = 100;
	let samps: i32 = 10;

	let cam: Ray = Ray{o:Vec3d{x:50.0,y:52.0,z:295.6}, d: Vec3d{x:0.0,y:-0.042612, z:-1.0}.normalise()};

	let cx: Vec3d = Vec3d{x:(width as f64)*0.5135/(height as f64),y:0.0,z:0.0};
	let cy: Vec3d = (cx%cam.d).normalise()*0.5135;
	let mut r: Vec3d = Vec3d{x:0.0,y:0.0,z:0.0};
	let mut c: Vec<Vec3d> = Vec::with_capacity((width*height) as usize);

	unsafe { c.set_len((width*height) as usize); }

	println!("Set c to len = {}", (width*height) as usize);
	println!("We set c.len = {}", c.len());

	for y in 0..height {
		println!("Rendering ({} spp) {}",samps*4,100*y/(height-1));
		
		for x in 0..width {
			let xi = y.pow(3);
			for sy in 0..2 {
				let i = ((height-y-1)*width+x) as usize;
				for sx in 0..2 {
					for s in 0..samps {
						let r1: f64 = 1.0; //2.0*rand::random::<f64>(); //erand48(xi);
						let r2: f64 = 1.0; //2.0*rand::random::<f64>(); //erand48(xi);

						let dx: f64;
						let dy: f64;

						if r1 < 1.0 {
							dx = r1.sqrt() - 1.0;
						} else {
							dx = 1.0 - (2.0-r1).sqrt();
						}

						if r2 < 1.0 {
							dy = r2.sqrt() - 1.0;
						} else {
							dy = 1.0 - (2.0-r2).sqrt();
						}

						let d: Vec3d = cx*((((sx as f64)+0.5 + dx)/2.0 + (x as f64))/ (width as f64) - 0.5) + cy*((((sy as f64)+0.5+dy)/2.0 + (y as f64))/(height as f64) - 0.5) + cam.d;
						println!("original dir  d.x: {}, d.y: {}, d.z: {}", d.x, d.y, d.z);
						let rad: Vec3d = radiance(&mut spheres, &Ray{o:cam.o+d*140.0,d:d.normalise()},0);
						//println!("rad.x: {}, rad.y: {}, rad.z: {}", rad.x, rad.y, rad.z);
						r = r + rad*(1.0/samps as f64);
						
					}
					let v: Vec3d = Vec3d{x: clamp(r.x), y: clamp(r.y), z: clamp(r.z)};
					//println!("v.x: {}, v.y: {}, v.z: {}", v.x, v.y, v.z);
					c[i] = c[i] + v*0.25;
					r = Vec3d{x:0.0,y:0.0,z:0.0};
				}
			}
		}
	}

	println!("Random = {}", rand::random::<f64>());
	println!("Random = {}", rand::random::<f64>());
	println!("Random = {}", rand::random::<f64>());



	// We create file options to write
	let file = OpenOptions::new().write(true).create(true).open("image.ppm").unwrap();

	// We create a buffered writer from the file we get
	let mut writer = BufWriter::new(&file);
	// Then we write to the file. write_all() calls flush() after the write as well.
	let mut b = format!("P3\n{} {}\n{}\n", width, height, 255).into_bytes();
	writer.write_all(&b);
	for i in 0..(width*height) as usize {
		b = format!("{} {} {}\n", to_int(c[i].x), to_int(c[i].y), to_int(c[i].z)).into_bytes();
		//b = format!("{} {} {}\n", c[i].x, c[i].y, c[i].z).into_bytes();
		//println!("c[i].x: {}, c[i].y: {}, c[i].z: {}", c[i].x, c[i].y, c[i].z);
		writer.write_all(&b);
	}

}



/*

 int main(int argc, char *argv[]){
   int w=1024, h=768, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
   Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
   Vec cx=Vec(w*.5135/h), cy=(cx%cam.d).norm()*.5135, r, *c=new Vec[w*h];
 #pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
   for (int y=0; y<h; y++){                       // Loop over image rows
     fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
     * 
     for (unsigned short x=0, Xi[3]={0,0,y*y*y}; x<w; x++)   // Loop cols
       for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
         for (int sx=0; sx<2; sx++, r=Vec()){        // 2x2 subpixel cols
           for (int s=0; s<samps; s++){
           * 
             double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
             double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
             Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                     cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
             r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps);
           } // Camera rays are pushed ^^^^^ forward to start in interior
           c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
         }
   }
   FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
   fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
   for (int i=0; i<w*h; i++)
     fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    */
