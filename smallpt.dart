#!/usr/local/dart/dart-sdk/bin/dart --checked

// smallpt.dart - a Dart port of 'smallpt' global illumination
// cf. http://www.kevinbeason.com/smallpt/

import 'dart:math';
import 'dart:io';

class Vec {
    double x, y, z;
    // Vec(this.x, this.y, this.z);
    Vec([double x = 0.0, double y = 0.0, double z = 0.0]) {
        this.x = x; this.y = y; this.z = z;
    }

    Vec operator+(Vec b) => new Vec(x+b.x, y+b.y, z+b.z);
    Vec operator-(Vec b) => new Vec(x-b.x, y-b.y, z-b.z);
    Vec operator*(double b) => new Vec(x*b, y*b, z*b);
    Vec mult(Vec b) => new Vec(x*b.x, y*b.y, z*b.z);
    Vec norm() {
        double l = sqrt(x * x + y * y + z * z);
        x /= l; y /= l; z /= l;
        return this;
    }
    double dot(Vec b) => x*b.x + y*b.y + z*b.z;
    Vec operator%(Vec b) =>
        new Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);
}

class Ray {
    Vec o, d;
    Ray(this.o, this.d);
}

// enum Refl_t
const int DIFF = 0;
const int SPEC = 1;
const int REFR = 2;

var rng = new Random();

class Sphere {
    double rad;
    Vec p, e, c;
    /* Refl_t */ int refl;
    Sphere(this.rad, this.p, this.e, this.c, this.refl);

    double intersect(Ray r) {
        Vec op = p - r.o;
        // double t;
        double eps = 1e-4;
        double b = op.dot(r.d);
        double det = b * b - op.dot(op) + rad * rad;
        if (det < 0) return /* 0 */ 1e20;
        det = sqrt(det);
        var t = b - det;
        if (t > eps) return t;
        t = b + det;
        if (t > eps) return t;
	return /* 0 */ 1e20;
    }
}

var spheres = new List<Sphere>();

double clamp(double x) => x<0 ? 0.0 : (x>1 ? 1.0 : x);
int toInt(double x) => (pow(x.clamp(0, 1), 1/2.2) * 255 + 0.5).toInt();

intersect(Ray r)
{
    Sphere id = null;
    const double inf = 1e20;
    double t = inf;
    for (var sphere in spheres) {
	double d = sphere.intersect(r);
	if (d < t) {
	    t = d;
	    id = sphere;
	}
    }
    return [t, id];
}

Vec radiance(Ray r, int depth)
{
    // var rng = new Random();
    double t;
    Sphere obj;

    var t_and_obj = intersect(r);
    t = t_and_obj[0];
    obj = t_and_obj[1];

    if (obj == null) return new Vec();
    Vec x = r.o + r.d * t;
    Vec n = (x - obj.p).norm();
    Vec nl = n.dot(r.d)<0 ? n : n*-1.0;
    Vec f = obj.c;
    double p = (f.x > f.y && f.x > f.z) ? f.x : (f.y > f.z ? f.y : f.z);
    depth++;

    if (depth > 5) {
	if (rng.nextDouble() < p)
	    f = f * (1/p);
	else return obj.e;
    }
    if (obj.refl == DIFF) {
	double r1 = 2 * PI * rng.nextDouble();
	double r2 = rng.nextDouble();
	double r2s = sqrt(r2);
	Vec w = nl;
	Vec u = ((w.x.abs()>0.1 ? new Vec(0.0, 1.0) : new Vec(1.0)) % w).norm();
	Vec v = w % u;
	Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
	return obj.e + f.mult(radiance(new Ray(x,d),depth));
    }
    if (obj.refl == SPEC)
	return obj.e + f.mult(radiance(new Ray(x,r.d-n*2.0*n.dot(r.d)),depth));

    var reflRay = new Ray(x, r.d - n*2.0*n.dot(r.d));
    bool into = n.dot(nl) > 0;
    double nc = 1.0, nt = 1.5;
    double nnt = into ? nc/nt : nt/nc;
    double ddn = r.d.dot(nl);
    double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
    if (cos2t < 0)
	return obj.e + f.mult(radiance(reflRay,depth));
    Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
    double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
    double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
    return obj.e + f.mult(depth > 2 ? (rng.nextDouble() < P ?
	radiance(reflRay, depth) * RP : radiance(new Ray(x,tdir),depth) * TP) :
	radiance(reflRay, depth) * Re + radiance(new Ray(x,tdir),depth) * Tr);
}

main()
{
    spheres.add(new Sphere(1e5,   new Vec( 1e5+1,40.8,81.6),   new Vec(),new Vec(.75,.25,.25),DIFF));//Left
    spheres.add(new Sphere(1e5,   new Vec(-1e5+99,40.8,81.6),  new Vec(),new Vec(.25,.25,.75),DIFF));//Rght
    spheres.add(new Sphere(1e5,   new Vec(50.0,40.8, 1e5),     new Vec(),new Vec(.75,.75,.75),DIFF));//Back
    spheres.add(new Sphere(1e5,   new Vec(50.0,40.8,-1e5+170), new Vec(),new Vec(),           DIFF));//Frnt
    spheres.add(new Sphere(1e5,   new Vec(50.0, 1e5, 81.6),    new Vec(),new Vec(.75,.75,.75),DIFF));//Botm
    spheres.add(new Sphere(1e5,   new Vec(50.0,-1e5+81.6,81.6),new Vec(),new Vec(.75,.75,.75),DIFF));//Top
    spheres.add(new Sphere(16.5,  new Vec(27.0,16.5,47.0),     new Vec(),new Vec(1.0,1.0,1.0)*.999, SPEC));//Mirr
    spheres.add(new Sphere(16.5,  new Vec(73.0,16.5,78.0),     new Vec(),new Vec(1.0,1.0,1.0)*.999, REFR));//Glas
    spheres.add(new Sphere(600.0, new Vec(50.0,681.6-.27,81.6),new Vec(12.0,12.0,12.0), new Vec(), DIFF)); //Lite

    List<String> args = (new Options()).arguments;
    int w = int.parse(args[0]);
    int h = int.parse(args[1]);
    int samps = (args.length >= 3) ? int.parse(args[2]) ~/ 4 : 1;
    Ray cam = new Ray(new Vec(50.0, 52.0, 295.6),
	(new Vec(0.0, -0.042612, -1.0)).norm());
    Vec cx = new Vec(w * 0.5135 / h);
    Vec cy = (cx % cam.d).norm() * 0.5135;
    // Vec r = null;

    var c = new List<Vec>(w * h);
    for (int i = 0; i < w * h; ++i)
	c[i] = new Vec();

    // var rng = new Random();

    for (int y = 0; y < h; ++y) {	// Loop over image rows
	stdout.write("\rRendering (${samps * 4} spp) ${(100.0 * y / (h-1)).toStringAsFixed(2)}%");
	for (int x = 0; x < w; ++x) {
	    int i = (h - y - 1) * w + x;
	    for	(int sy = 0; sy < 2; ++sy) {
		Vec r = new Vec();
		for (int sx = 0; sx < 2; ++sx) {
		    for (int s = 0; s < samps; ++s) {
			double r1 = 2 * rng.nextDouble();
			double dx = (r1<1) ? sqrt(r1) - 1 : 1 - sqrt(2-r1);
			double r2 = 2 * rng.nextDouble();
			double dy = (r2<1) ? sqrt(r2) - 1 : 1 - sqrt(2-r2);
			Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
			    cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
			r = r +
			    radiance(new Ray(cam.o + d*140.0, d.norm()),0) *
			    (1 / samps);
		    }
		    c[i] = c[i] +
			new Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
		}
	    }
	}
    }

    var image_file = new File("image.ppm");
    var sink = image_file.openWrite();
    sink.write("P3\n${w} ${h}\n255\n");
    for (var v in c) {
	int r = toInt(v.x);
	int g = toInt(v.y);
	int b = toInt(v.z);
	sink.write("${r} ${g} ${b} ");
    }
    sink.close();

    print("\nend of main");
}

// eof	vim:set syntax=dart:
