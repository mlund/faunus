#include "povray.h"

povray::povray() {
  anionTxt   = "greyish";
  cationTxt  = "redish"; 
  neutralTxt = "white";
  defineTextures();
};

void povray::defineTextures() {
  ostr << \
    "#declare white=texture {\n"\
    " pigment {color rgb <1,1,1>}\n"\
    " finish {phong .9 ambient .1 reflection 0.2}\n"\
    "}\n"\
    "#declare black=texture {\n"\
    " pigment {color rgb <0,0,0>}\n"\
    " finish {phong .9 ambient .2 reflection .2}\n"\
    "}\n"\
    "#declare transp=texture {\n"\
    " pigment {color rgbf <1,1,1,.9>}\n"\
    " finish {phong .1 ambient .2}\n"\
    "}\n"\
    "#declare greyish=texture {\n"\
    " pigment {color rgbf <.5,.5,.5,.7>}\n"\
    " finish {phong .9 ambient .1 reflection .2}\n"\
    "}\n" \
    "#declare redish=texture {\n"\
    " pigment {color rgb <1,0,0>}\n"\
    " finish {phong .9 ambient .1 reflection .2}\n"\
    "}\n" \
    "#declare light_object = sphere { 0,2 pigment { rgbf <1,1,1,1> } finish { specular 5 metallic } }\n";
};

void povray::coordinates(float r) {
  ostr << \
    "union {\n"\
    " cylinder {\n  <0,0,"<<-r<<">,<0,0,"<<r<<">,.5\n"\
    "  pigment {color rgb <0,.6,.3>}\n"\
    "  finish {ambient .3}\n"\
    " }\n"\
    " cylinder {\n  <0,"<<-r<<",0>,<0,"<<r<<",0>,.5\n"\
    "  pigment {color rgb <0,.6,.3>}\n"\
    "  finish {ambient .3}\n"\
    " }\n"\
    " cylinder {\n  <"<<-r<<",0,0>,<"<<r<<",0,0>,.5\n"\
    "  pigment {color rgb <0,.6,.3>}\n"\
    "  finish {ambient .3}\n"\
    " }\n"\
    "}\n";
};

void povray::zaxis(double r) {
  ostr << "cylinder {\n"
       << "  <0,0,"<<-r<<">, <0,0,"<<r<<">, 0.5\n"
       << "  texture {transp}\n"
       << "}\n";
};

void povray::glowingsphere(particle &p) {
  ostr << "light_source {<"
    << p.x << "," << p.y << "," << p.z
    << "> color rgb <1,0.75,0> looks_like { light_object }}\n";
}

void povray::sphere(double r) {
  //camera position
  ostr << "camera {" << endl
       << "  location <100,100,0>" << endl
       << "  look_at <0,0,0>\n}" << endl;

  //light sources
  ostr << "light_source {\n"
       << "  <"<<r*1.5<<","<<r*1.5<<","<<r*0.5<<">\n"
       << "  color rgb <1,1,1>\n}\n";
  ostr << "light_source {\n"
       << "  <"<<-r*1.5<<","<<r*1.5<<","<<-r*0.5<<">\n"
       << "  color rgb <1,1,1>\n}" << endl;

  ostr << "sphere {<0,0,0>," << r
       << " texture {transp}"
       << "}\n";
};

void povray::bond(point &p1, point &p2, float radius) {
  ostr << " cylinder {\n"
       << "  <"<<p1.x<<","<<p1.y<<","<<p1.z<<">,"
       << "  <"<<p2.x<<","<<p2.y<<","<<p2.z<<">,"<<radius
       << "  texture {"<<neutralTxt<<"}\n }\n";
};

void povray::add(vector<particle> &p, group &g) {
  string txt;
  ostr << setiosflags( ios::fixed );
  ostr.precision(1);
  ostr << "union {\n";
  for (int i=g.beg; i<=g.end; i++) {
    if (p[i].charge>0)  txt=cationTxt;
    if (p[i].charge<0)  txt=anionTxt;
    if (p[i].charge==0) txt=neutralTxt;

    if (p[i].charge!=0 && p[i].radius==1.8 && 1>2) {
      glowingsphere(p[i]);
    } else if (p[i].radius>0) {
      ostr << " sphere {<"<<p[i].x<<","<<p[i].y<<","<<p[i].z<<">,"<<p[i].radius
	   << " texture {"<<txt<<"}}" << endl;
    };

    //chain monomers...
    if (g.chain==true) {
      if (i==g.beg && g.graftpoint>-1)
        bond( p[i], p[g.graftpoint] );
      if (i!=g.end)
        bond( p[i], p[i+1] );
    };
  };
  ostr << "}\n";
};

void povray::cube(double boxlen) {
  double a=boxlen/2.;
  point p1,p2;
  ostr << "union {\n";

  p1.x=a;  p1.y=a; p1.z=a;
  p2.x=-a; p2.y=a; p2.z=a;
  cylinder(p1,p2);
  p1.x=a; p1.y=a; p1.z=a;
  p2.x=a; p2.y=-a; p2.z=a;
  cylinder(p1,p2);
  p1.x=a; p1.y=a; p1.z=a;
  p2.x=a; p2.y=a; p2.z=-a;
  cylinder(p1,p2);
  p1.x=-a; p1.y=-a; p1.z=-a;
  p2.x=a; p2.y=-a; p2.z=-a;
  cylinder(p1,p2);
  p1.x=-a; p1.y=-a; p1.z=-a;
  p2.x=-a; p2.y=a; p2.z=-a;
  cylinder(p1,p2);
  p1.x=-a; p1.y=-a; p1.z=-a;
  p2.x=-a; p2.y=-a; p2.z=a;
  cylinder(p1,p2);  
  ostr << "}\n";
  p1.x=a; p1.y=-a; p1.z=a;
  p2.x=-a; p2.y=-a; p2.z=a;
  cylinder(p1,p2);
  p1.x=a; p1.y=-a; p1.z=a;
  p2.x=a; p2.y=-a; p2.z=-a;
  cylinder(p1,p2);
  p1.x=-a; p1.y=a; p1.z=a;
  p2.x=-a; p2.y=a; p2.z=-a;
  cylinder(p1,p2);
  p1.x=-a; p1.y=a; p1.z=a;
  p2.x=-a; p2.y=-a; p2.z=a;
  cylinder(p1,p2);
  p1.x=a; p1.y=a; p1.z=-a;
  p2.x=-a; p2.y=a; p2.z=-a;
  cylinder(p1,p2);
  p1.x=a; p1.y=a; p1.z=-a;
  p2.x=a; p2.y=-a; p2.z=-a;
  cylinder(p1,p2);
};

void povray::cylinder(point &p1, point &p2, double r) {
  ostr << 
    " cylinder { <" 
       << p1.x <<","<<p1.y<<","<<p1.z<<">,<"
       << p2.x <<","<<p2.y<<","<<p2.z<<">,"<<r<<"\n"\
    "  pigment {color rgb <0,.6,.3>}\n"\
    "  finish {ambient .3}\n"\
    " }\n";
};


//Write POVRAY data to file
void povray::save(string filename) {
  ofstream f( filename.c_str() );
  if (f) {
    f << ostr.str();
    f.close();
    //cout << "*** POVRAY data written to '"<<filename<<"'"<< endl;
  } else cout << "*** Error creating POVRAY file." << endl;
};
