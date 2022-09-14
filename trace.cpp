#include "trace.H"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <getopt.h>
#ifdef __APPLE__
#define MAX std::numeric_limits<double>::max()
#else
#include <values.h>
#define MAX DBL_MAX
#endif

using namespace std;
// return the determinant of the matrix with columns a, b, c.
double det(const SlVector3 &a, const SlVector3 &b, const SlVector3 &c) {
    return a[0] * (b[1] * c[2] - c[1] * b[2]) +
        b[0] * (c[1] * a[2] - a[1] * c[2]) +
            c[0] * (a[1] * b[2] - b[1] * a[2]);
}

inline double sqr(double x) {return x*x;}

bool Triangle::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    SlVector3 ba = a-b;
    SlVector3 ca = a-c;
    SlVector3 ea = a-r.e;
    double detA = det(ba, ca, r.d);
    double t = det(ba, ca, ea)/detA;
    if (t < t0 || t > t1) return false;
    double beta = det(ea, ca, r.d)/detA;
    if (beta < 0 || beta > 1) return false;
    double gamma = det(ba, ea, r.d)/detA;
    if (gamma < 0.0 || gamma > 1.0-beta) return false;
    hr.t = t;
    hr.p = r.e + t * r.d;
    hr.n = cross(ba,ca);
    normalize(hr.n);
    hr.alpha = 1.0 - beta - gamma;
    hr.beta = beta;
    hr.gamma = gamma;
    return true;
}

bool TrianglePatch::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    bool temp = Triangle::intersect(r,t0,t1,hr);
    if (temp) {
        hr.n = hr.alpha * n1 + hr.beta * n2 + hr.gamma * n3;
        normalize(hr.n);
    }
    return temp;
}


bool Sphere::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    double ddotemc = dot(r.d, r.e-c);
    double d2 = sqrMag(r.d);

    double disc = sqr(ddotemc) - d2 * (sqrMag(r.e-c) - rad*rad);

    if (disc < 0) return false;
    double root1 = (-ddotemc + sqrt(disc)) / d2;
    double root2 = (-ddotemc - sqrt(disc)) / d2;

    double t = root1;
    if (root1 < 0 || (root2 > 0 && root2 < root1)) t = root2;
    if (t < t0 || t > t1) return false;
  
    hr.t = t;
    hr.p = r.e + t * r.d;
    hr.n = (hr.p - c) / rad;
    return true;
}


Tracer::Tracer(const std::string &fname) {
    std::ifstream in(fname.c_str(), std::ios_base::in);
    std::string line;
    char ch;
    Fill fill;
    bool coloredlights = false;
    while (in) {
        getline(in, line);
        switch (line[0]) {
            case 'b': {
                std::stringstream ss(line);
                ss>>ch>>bcolor[0]>>bcolor[1]>>bcolor[2];
                break;
            }

            case 'v': {
                getline(in, line);
                std::string junk;
                std::stringstream fromss(line);
                fromss>>junk>>eye[0]>>eye[1]>>eye[2];

                getline(in, line);
                std::stringstream atss(line);
                atss>>junk>>at[0]>>at[1]>>at[2];

                getline(in, line);
                std::stringstream upss(line);
                upss>>junk>>up[0]>>up[1]>>up[2];

                getline(in, line);
                std::stringstream angless(line);
                angless>>junk>>angle;

                getline(in, line);
                std::stringstream hitherss(line);
                hitherss>>junk>>hither;

                getline(in, line);
                std::stringstream resolutionss(line);
                resolutionss>>junk>>res[0]>>res[1];
                break;
            }

            case 'p': {
                bool patch = false;
                std::stringstream ssn(line);
                unsigned int nverts;
                if (line[1] == 'p') {
                    patch = true;
                    ssn>>ch;
                }
                ssn>>ch>>nverts;
                std::vector<SlVector3> vertices;
                std::vector<SlVector3> normals;
                for (unsigned int i=0; i<nverts; i++) {
                    getline(in, line);
                    std::stringstream ss(line);
                    SlVector3 v,n;
                    if (patch) ss>>v[0]>>v[1]>>v[2]>>n[0]>>n[1]>>n[2];
                    else ss>>v[0]>>v[1]>>v[2];
                    vertices.push_back(v);
                    normals.push_back(n);
                }
                bool makeTriangles = false;
                if (vertices.size() == 3) {
                    if (patch) {
                        surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
                        normals [0], normals [1], normals [2]), fill));
                    } else {
                        surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                    }
                } else if (vertices.size() == 4) {
                    SlVector3 n0 = cross(vertices[1] - vertices[0], vertices[2] - vertices[0]);
                    SlVector3 n1 = cross(vertices[2] - vertices[1], vertices[3] - vertices[1]);
                    SlVector3 n2 = cross(vertices[3] - vertices[2], vertices[0] - vertices[2]);
                    SlVector3 n3 = cross(vertices[0] - vertices[3], vertices[1] - vertices[3]);
                    if (dot(n0,n1) > 0 && dot(n0,n2) > 0 && dot(n0,n3) > 0) {
                        makeTriangles = true;
                        if (patch) {
                            surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
                            normals[0], normals[1], normals[2]), fill));
                            surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[2], vertices[3],
                            normals[0], normals[2], normals[3]), fill));
                        } else {
                            surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                            surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[2], vertices[3]), fill));
                        }
                    }
                    if (!makeTriangles) {
                        std::cerr << "I didn't make triangles.  Poly not flat or more than quad.\n";
                    }
                }
                break;
            }

            case 's' : {
                std::stringstream ss(line);
                SlVector3 c;
                double r;
                ss>>ch>>c[0]>>c[1]>>c[2]>>r;
                surfaces.push_back(std::pair<Surface *, Fill>(new Sphere(c,r), fill));
                break;
            }
      
            case 'f' : {
                std::stringstream ss(line);
                ss>>ch>>fill.color[0]>>fill.color[1]>>fill.color[2]>>fill.kd>>fill.ks>>fill.shine>>fill.t>>fill.ior;
                break;
            }

            case 'l' : {
                std::stringstream ss(line);
                Light l;
                ss>>ch>>l.p[0]>>l.p[1]>>l.p[2];
                if (!ss.eof()) {
                    ss>>l.c[0]>>l.c[1]>>l.c[2];
                    coloredlights = true;
                }
                lights.push_back(l);
                break;
            }

            default:
            break;
        }
    }
    if (!coloredlights) for (unsigned int i=0; i<lights.size(); i++) lights[i].c = 1.0/sqrt(lights.size());
    im = new SlVector3[res[0]*res[1]];
    shadowbias = 1e-6;
    samples = 1;
    aperture = 0.0;
}

Tracer::~Tracer() {
    if (im) delete [] im;
    for (unsigned int i=0; i<surfaces.size(); i++) delete surfaces[i].first;
}

inline SlVector3 reflect(const SlVector3 &v, const SlVector3 &n) {
    return 2*dot(v,n)*n-v;
};

inline bool refract(const SlVector3 &v, const SlVector3 &n, double ir, SlVector3 &r) {
    double tmp = 1.0 - (1.0 - sqr(dot(v,n))) / sqr(ir);
    if (tmp < 0.0) return false;
    r =  (n * dot(v,n) - v) / ir  -  sqrt(tmp) * n;
    return true;
}

SlVector3 Tracer::shade(const HitRecord &hr) const {
    if (color) return hr.f.color;

    SlVector3 color(0.0);
    HitRecord dummy;
    HitRecord rhr;
    bool reflect=false;
    float t1=10000;
    SlVector3 viewer=hr.v;
    SlVector3 position=hr.p;
    normalize(viewer);
    SlVector3 normal2=hr.n;
    SlVector3 reflectionE=2*normal2*dot(viewer,normal2)-hr.v;
    normalize(reflectionE);
    
    Ray reflectionRayE(hr.p,reflectionE);

    
    
    for (unsigned int i=0; i<lights.size(); i++) {
        const Light &light = lights[i];
        bool shadow = false;
        
        SlVector3 position=hr.p;
        SlVector3 LightPos=light.p;
        SlVector3 light_po2=LightPos-position;
        normalize(light_po2);
        
        SlVector3 reflection2=2*normal2*dot(light_po2,normal2)-light_po2;
        normalize(reflection2);
        Ray reflectionRay(hr.p,reflection2);
        
        
       
        
        // Step 3 Check for shadows here
        SlVector3 surfacetoLight=light.p-hr.p;
        normalize(surfacetoLight);
        Ray surfacetoLightRay = Ray(hr.p,surfacetoLight);
        HitRecord newhr;
        
        for(int j=0;j<surfaces.size();j++){
            if(surfaces[j].first->intersect(surfacetoLightRay,0.00001,100000,newhr)){
                shadow=true;
            }
        }

        if (!shadow) {
 
            // Step 2 do shading here
            SlVector3 position=hr.p;
            SlVector3 LightPos=light.p;
            SlVector3 normal=hr.n;
            normalize(normal);
            SlVector3 light_po=LightPos-position;
            normalize(light_po);
            float dfd=dot(light_po,normal);
            if(dfd<0){
                dfd=0;
            }
            double diffuse=hr.f.kd*dfd;
            SlVector3 reflection=2*normal*(dot(light_po,normal))-light_po;
            normalize(reflection);
            float ddd=dot(reflection,viewer);
            if(ddd<0){
                ddd=0;
            }
            double specular=hr.f.ks*pow(ddd,hr.f.shine);
            color+=(diffuse+specular+0.1)*light.c*hr.f.color;
            
            
        }
    }
    
    
    
    // Step 4 Add code for computing reflection color here
    for(int j=0;j<surfaces.size();j++){
        if(surfaces[j].first->intersect(reflectionRayE,0.0001,t1,rhr)){
            rhr.f = hr.f;
            rhr.raydepth = hr.raydepth+1;
            rhr.v = reflectionRayE.e - hr.p;
            normalize(rhr.v);
            normalize(rhr.n);
            cout<<reflect<<endl;
            reflect=true;
        }
    }
    SlVector3 reflectColor;
    if(rhr.raydepth<5){
        reflectColor=shade(rhr);
        cout<<"111   "<<rhr.raydepth<<endl;
    }
    color+=reflectColor*hr.f.ks;
    // Step 5 Add code for computing refraction color here
    //color = hr.f.color;
    return color;
}

SlVector3 Tracer::trace(const Ray &r, double t0, double t1) const {
    HitRecord hr;
    SlVector3 color(bcolor);
  
    bool hit = false;
    for (unsigned int k=0; k<surfaces.size(); k++) {
        const std::pair<Surface *, Fill> &s  = surfaces[k];
        if (s.first->intersect(r, t0, t1, hr)) {
            t1 = hr.t;
            hr.f = s.second;
            hr.raydepth = r.depth;
            hr.v = r.e - hr.p;
            normalize(hr.v);
            normalize(hr.n);
            hit = true;
        }
    }

    if (hit) color = shade(hr);
    return color;
}

void Tracer::traceImage() {
    // set up coordinate system
    SlVector3 w = eye - at;
    w /= mag(w);
    SlVector3 u = cross(up,w);
    normalize(u);
    SlVector3 v = cross(w,u);
    normalize(v);

    double d = mag(eye - at);
    double h = tan((M_PI/180.0) * (angle/2.0)) * d;
    double l = -h;
    double r = h;
    double b = h;
    double t = -h;

    SlVector3 *pixel = im;

    for (unsigned int j=0; j<res[1]; j++) {
        for (unsigned int i=0; i<res[0]; i++, pixel++) {

            SlVector3 result(0.0,0.0,0.0);

            for (int k = 0; k < samples; k++) {

                double rx = 1.1 * rand() / RAND_MAX;
                double ry = 1.1 * rand() / RAND_MAX;

                double x = l + (r-l)*(i+rx)/res[0];
                double y = b + (t-b)*(j+ry)/res[1];
                SlVector3 dir = -d * w + x * u + y * v;
    
                Ray r(eye, dir);
                normalize(r.d);

                result += trace(r, hither, MAX);

            }
            (*pixel) = result / samples;
        }
    }
}

void Tracer::writeImage(const std::string &fname) {
#ifdef __APPLE__
    std::ofstream out(fname, std::ios::out | std::ios::binary);
#else
    std::ofstream out(fname.c_str(), std::ios_base::binary);
#endif
    out<<"P6"<<"\n"<<res[0]<<" "<<res[1]<<"\n"<<255<<"\n";
    SlVector3 *pixel = im;
    char val;
    for (unsigned int i=0; i<res[0]*res[1]; i++, pixel++) {
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[0])) * 255.0);
        out.write (&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[1])) * 255.0);
        out.write (&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[2])) * 255.0);
        out.write (&val, sizeof(unsigned char));
    }
    out.close();
}


int main(int argc, char *argv[]) {
    int c;
    double aperture = 0.0;
    int samples = 1;
    int maxraydepth = 5;
    bool color = false;
    bool casting = false;
    while ((c = getopt(argc, argv, "a:s:d:cb")) != -1) {
        switch(c) {
            case 'a':
            aperture = atof(optarg);
            break;
            case 's':
            samples = atoi(optarg);
            break;
            case 'c':
            color = true;
            break;
            case 'b':
            casting = true;
            break;
            case 'd':
            maxraydepth = atoi(optarg);
            break;
            default:
            abort();
        }
    }
    /*
    if (argc-optind != 2) {
        std::cout<<"usage: trace [opts] input.nff output.ppm"<<std::endl;
        for (unsigned int i=0; i<argc; i++) std::cout<<argv[i]<<std::endl;
        //exit(0);
    }
     */
    //Tracer tracer(argv[optind++]);
    Tracer tracer("/Users/billlu/Desktop/Study/Others/Graphic/as2-StarterVersion/InputFiles/mount.nnf");
    tracer.aperture = aperture;
    tracer.samples = samples;
    tracer.color = color;
    tracer.maxraydepth = maxraydepth;
    tracer.traceImage();
    //tracer.writeImage(argv[optind++]);
    tracer.writeImage("/Users/billlu/Desktop/Study/Others/Graphic/as2-StarterVersion/InputFiles/mount.ppm");
};
