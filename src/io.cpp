#include "io.h"


bool Faunus::IO::readFile(const std::string &file, std::vector<std::string> &v) {
    std::ifstream f( file );
    if (f) {
        std::string s;
        while (getline(f,s))
            v.push_back(s);
        f.close();
        return true;
    }
    std::cerr << "# WARNING! FILE " << file << " NOT READ!\n";
    return false;
}

bool Faunus::IO::writeFile(const std::string &file, const std::string &s, std::ios_base::openmode mode) {
    std::ofstream f(file, mode);
    if (f) {
        f << s;
        return true;
    }
    return false;
}

void Faunus::IO::strip(std::vector<std::string> &v, const std::string &pat) {
    auto iter=v.begin();
    while (iter!=v.end())
        if ((*iter).find(pat) != std::string::npos)
            v.erase(iter);
        else ++iter;
}

int Faunus::FormatXTC::getNumAtoms() { return natoms_xtc; }

bool Faunus::FormatAAM::keepcharges = true;

bool Faunus::FormatXTC::open(std::string s) {
    if (xd!=NULL)
        close();
    xd = xdrfile_open(&s[0], "r");
    if (xd!=NULL) {
        int rc = read_xtc_natoms(&s[0], &natoms_xtc); // get number of atoms
        if (rc==exdrOK) {
            x_xtc = new rvec [natoms_xtc]; // resize coordinate array
            return true;
        }
    } else
        std::cerr << "# ioxtc error: xtc file could not be opened." << endl;
    return false;
}

void Faunus::FormatXTC::close() {
    xdrfile_close(xd);
    xd=NULL;
    delete[] x_xtc;
}

Faunus::FormatXTC::FormatXTC(double len) {
    prec_xtc = 1000.;
    time_xtc=step_xtc=0;
    setbox(len);
    xd=NULL;
    x_xtc=NULL;
}

Faunus::FormatXTC::~FormatXTC() {
    close();
}

void Faunus::FormatXTC::setbox(double x, double y, double z) {
    assert(x>0 && y>0 && z>0);
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            xdbox[i][j]=0;
    xdbox[0][0]=0.1*x; // corners of the
    xdbox[1][1]=0.1*y; // rectangular box
    xdbox[2][2]=0.1*z; // in nanometers!
}

void Faunus::FormatXTC::setbox(double len) { setbox(len,len,len); }

void Faunus::FormatXTC::setbox(const Faunus::Point &p) { setbox(p.x(), p.y(), p.z()); }

std::string Faunus::FormatPQR::writeCryst1(const Faunus::Point &len, const Faunus::Point &angle) {
    char buf[500];
    sprintf(buf, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
            len.x(),len.y(),len.z(),angle.x(),angle.y(),angle.z());
    return std::string(buf);
}
