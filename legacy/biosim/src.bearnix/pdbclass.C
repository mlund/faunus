#include "pdbclass.h"

pdbclass::pdbclass() {
  modelCnt=1;
  remarkCnt=1;
};

void pdbclass::remark(double s) {
  out.setf(ios::left, ios::adjustfield);
  out << "REMARK "; //1-7
  out << setw(3) << remarkCnt++ << " "; //8-10+1
  out << s << endl; //12-
};

void pdbclass::conect(int a, int b) {
  out << "CONECT";        //1-6
  out.setf(ios::right, ios::adjustfield);
  out << setw(5) << a;    //7-11
  out << setw(5) << b;
  out << endl;
};

void pdbclass::hetatm(format &f) {
  out.setf(ios::left, ios::adjustfield);
  out << "ATOM  ";                             // 1-6
  out.setf(ios::right, ios::adjustfield);
  out << setw(5) << f.atomnumber << " ";       // 7-11 + 1
  out.setf(ios::left, ios::adjustfield);
  out << setw(4) << f.atomname << " ";         // 13-16 + 1
  out << setw(3) << f.residuename << " ";      // 18-20 + 1
  out << setw(1) << f.chainid;                 // 22
  out.setf(ios::right, ios::adjustfield);
  out << setw(4) << f.residuenumber << "    "; // 23-26 + 4
  out.precision(3);
  out.setf(ios::fixed, ios::floatfield);
  out << setw(8) << f.x;  //31-38
  out << setw(8) << f.y;  //39-46
  out << setw(8) << f.z;  //47-54
  out << "                  "; //-72
  out.setf(ios::left, ios::adjustfield);
  out << setw(4) << f.seqid;
  out << endl;
};

void pdbclass::atom(format &f) {
  //out << "12345678901234567890123456789012345678901234567890" << endl;
  out.setf(ios::left, ios::adjustfield);
  out << "ATOM  ";                             // 1-6
  out.setf(ios::right, ios::adjustfield);
  out << setw(5) << f.atomnumber << " ";       // 7-11 + 1
  out.setf(ios::left, ios::adjustfield);
  out << setw(4) << f.atomname << " ";         // 13-16 + 1
  out << setw(3) << f.residuename << " ";      // 18-20 + 1
  out << setw(1) << f.chainid;                 // 22
  out.setf(ios::right, ios::adjustfield);
  out << setw(4) << f.residuenumber << "    "; // 23-26 + 4
  out.precision(3);
  out.setf(ios::fixed, ios::floatfield);
  out << setw(8) << f.x;  //31-38
  out << setw(8) << f.y;  //39-46
  out << setw(8) << f.z;  //47-54
//PQR
  out << setw(8) << " " << f.charge << " " << f.radius << " ";
//PQR
  out << "                  "; //-72
  out.setf(ios::left, ios::adjustfield);
  out << setw(4) << f.seqid;
  out << endl;
};

void pdbclass::model(int i) {
  if (i==0) i=modelCnt++;
  out.setf(ios::right, ios::adjustfield);
  out << "MODEL" << "    ";     // 1-6 + 4
  out << setw(4) << i << endl; // 11-14
};

void pdbclass::endmdl() {
  out << "ENDMDL" << endl;
};

void pdbclass::save(string filename) {
  ofstream f( filename.c_str() );
  if (f) {
    f << out.str();
    f.close();
  } else cout << "pdb file error!!\n";
};

void pdbclass::load_particles(vector<particle> &p) {
  format f;
  f.atomnumber=0;
  f.atomname="O";
  f.residuename="UNK";
  f.chainid=" ";
  f.residuenumber=0;
  model();
  for (int i=0; i<p.size(); i++) {
    f.atomnumber=i+1;
    f.residuenumber=i+1;
    f.x = p[i].x;
    f.y = p[i].y;
    f.z = p[i].z;
    f.charge=0;
    if (p[i].charge>0) f.charge++;
    if (p[i].charge<0) f.charge--;
    f.radius = p[i].radius;
    atom(f);
  };
  endmdl();
};
