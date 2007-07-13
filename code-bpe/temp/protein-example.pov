#declare white=texture {
 pigment {color rgb <1,1,1>}
 finish {phong .9 ambient .1 reflection 0.2}
}
#declare black=texture {
 pigment {color rgb <0,0,0>}
 finish {phong .9 ambient .2 reflection .2}
}
#declare transp=texture {
 pigment {color rgbf <1,1,1,.9>}
 finish {phong .1 ambient .2}
}
#declare greyish=texture {
 pigment {color rgbf <.5,.5,.5,.7>}
 finish {phong .9 ambient .1 reflection .2}
}
#declare redish=texture {
 pigment {color rgb <1,0,0>}
 finish {phong .9 ambient .1 reflection .2}
}
light_source {
  <150,150,50>
  color rgb <1,1,1>
}
light_source {
  <-150,150,-50>
  color rgb <1,1,1>
}
camera {
  location <100,100,0>
  look_at <0,0,0>
}
cylinder {<0,0,0>,<0,0,100>,50 texture {transp}}
cylinder {<0,0,0>,<0,0,100>,0.5 texture {transp}}
sphere {<-11.827,-2.37029,7.99883>,2 texture {white}}
sphere {<-10.897,-0.937295,8.22783>,3.7 texture {white}}
sphere {<-8.66704,-1.88729,3.93783>,3.3 texture {white}}
sphere {<-5.52704,-3.38729,3.96783>,3.4 texture {white}}
sphere {<-3.88704,0.952705,1.72783>,3.8 texture {white}}
sphere {<-8.40704,2.09271,5.71783>,3.8 texture {white}}
sphere {<-5.23704,-0.247295,8.42783>,3.6 texture {white}}
sphere {<-0.107044,1.15271,5.62783>,3.7 texture {white}}
sphere {<-2.92704,5.05271,6.99783>,2.9 texture {white}}
sphere {<-4.61704,4.53271,10.3078>,3.6 texture {white}}
sphere {<1.11296,2.66271,11.0178>,3.9 texture {white}}
sphere {<1.22296,6.15271,7.17783>,3.8 texture {white}}
sphere {<-2.18704,9.58271,10.6378>,3.7 texture {white}}
sphere {<-1.42704,8.03271,14.5278>,4.1 texture {white}}
sphere {<4.48296,6.92271,12.8778>,3.1 texture {white}}
sphere {<4.99296,10.4527,11.8878>,3.1 texture {white}}
sphere {<2.79296,12.2327,15.2778>,3.7 texture {white}}
sphere {<7.51296,10.4527,17.7278>,3.8 texture {white}}
sphere {<10.813,9.81271,14.8778>,2.9 texture {white}}
sphere {<10.923,6.08271,12.5178>,3.6 texture {white}}
sphere {<8.58296,8.27271,10.1378>,3.4 texture {white}}
sphere {<8.42296,4.51271,8.96783>,3.6 texture {white}}
sphere {<8.79296,2.54271,13.0178>,3.8 texture {white}}
sphere {<4.69296,2.80271,16.2078>,3.6 texture {white}}
sphere {<6.01296,4.95271,19.6378>,3.3 texture {white}}
sphere {<5.13296,3.11271,24.4978>,3.7 texture {white}}
sphere {<3.36296,8.13271,23.3678>,3.8 texture {white}}
sphere {<3.37296,7.08271,17.8878>,3.8 texture {white}}
sphere {<0.462956,2.32271,19.9078>,3.6 texture {white}}
sphere {<-1.66704,5.68271,24.0278>,3.7 texture {white}}
sphere {<-2.14704,8.45271,20.1978>,3.6 texture {white}}
sphere {<-3.53704,4.47271,17.1278>,3.6 texture {white}}
sphere {<-4.88704,2.13271,20.2678>,3.6 texture {white}}
sphere {<-6.49704,6.39271,23.0978>,3.8 texture {white}}
sphere {<-8.26704,7.89271,18.7278>,3.5 texture {white}}
sphere {<-8.68704,5.36271,15.5878>,3.8 texture {white}}
sphere {<-8.44704,0.902705,16.9978>,3.9 texture {white}}
sphere {<-10.607,2.20271,22.0578>,3.4 texture {white}}
sphere {<-10.907,-2.22729,22.8578>,3.3 texture {white}}
sphere {<-7.27704,-3.25729,21.4178>,3.6 texture {white}}
sphere {<-6.05704,0.122705,25.2078>,3.6 texture {white}}
sphere {<-9.41704,-0.397295,28.2278>,3.7 texture {white}}
sphere {<-7.84704,-5.17729,28.3278>,2.9 texture {white}}
sphere {<-5.99704,-8.16729,29.0378>,3.4 texture {white}}
sphere {<-3.09704,-6.23729,28.0278>,3.3 texture {white}}
sphere {<-0.807044,-3.28729,29.5678>,3.5 texture {white}}
sphere {<0.282956,0.142705,26.2278>,3.6 texture {white}}
sphere {<3.65296,-1.10729,28.7478>,3.6 texture {white}}
sphere {<1.37296,-6.08729,28.7178>,3.8 texture {white}}
sphere {<0.302956,-4.59729,23.6478>,3.6 texture {white}}
sphere {<5.57296,-1.79729,22.5378>,3.9 texture {white}}
sphere {<7.37296,-5.94729,25.2778>,3.8 texture {white}}
sphere {<2.75296,-9.03729,21.9178>,3.8 texture {white}}
sphere {<3.69296,-5.70729,19.0678>,3.6 texture {white}}
sphere {<9.21296,-4.80729,19.7978>,3.6 texture {white}}
sphere {<8.41296,-8.61729,17.3278>,3.7 texture {white}}
sphere {<13.263,-6.30729,16.9778>,3.6 texture {white}}
sphere {<13.813,-4.19729,20.4078>,2.9 texture {white}}
sphere {<13.613,-1.59729,17.6278>,3.6 texture {white}}
sphere {<11.213,0.582705,20.0078>,2.9 texture {white}}
sphere {<10.433,1.55271,16.9278>,3.8 texture {white}}
sphere {<6.03296,-0.797295,15.9778>,3.4 texture {white}}
sphere {<6.53296,-2.15729,12.2378>,3.3 texture {white}}
sphere {<2.14296,-2.27729,8.80783>,3.9 texture {white}}
sphere {<4.93296,-6.37729,8.67783>,3.8 texture {white}}
sphere {<6.12296,-5.39729,14.4078>,3.8 texture {white}}
sphere {<-0.107044,-2.88729,13.6478>,3.9 texture {white}}
sphere {<0.402956,-6.99729,10.2178>,3.8 texture {white}}
sphere {<-0.587044,-9.10729,15.3178>,3.4 texture {white}}
sphere {<-3.24704,-5.73729,16.6178>,3.6 texture {white}}
sphere {<-5.15704,-6.43729,12.3378>,3.4 texture {white}}
sphere {<-7.02704,-11.1073,14.6178>,3.7 texture {white}}
sphere {<-7.24704,-8.13729,17.9478>,3.7 texture {white}}
sphere {<-8.44704,-4.27729,14.3778>,3.6 texture {white}}
sphere {<-10.637,-6.88729,11.7678>,3.3 texture {white}}
sphere {<-12.977,-9.39729,14.9378>,3.9 texture {white}}
