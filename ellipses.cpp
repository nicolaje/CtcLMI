#include <ibex/ibex.h>

#include "vibes.h"


#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "CtcLMI.h"
#include "CtcBoxLMI.h"
#include "vibes/vibes.h"
#include "CtcEllipsoid.h"
#include "CtcEllipsoids.h"

#define MAX_BISSECT 100

using namespace ibex;

void contractAndStore(IntervalVector &box, Ctc &contractor, vector<IntervalVector> &storage);
void exemple();
void invertEllipsoid();
void invertLMI();
void triangleBall();
void ctcAndDraw1();
void ctcAndDraw2();
void ctcAndDraw3();
void sivia(IntervalVector search, Ctc &cIn, Ctc &cOut, vector<IntervalVector> &in, vector<IntervalVector> &out, vector<IntervalVector> &perhaps);

int main()
{
    vibes::beginDrawing();
    //    exemple();
    NumConstraint cNu("x1", "x2", "x3", "(x1+x2)*x3>=1");
    NumConstraint cNu2("x1", "x2", "x3", "(x1+x2)*x3<=3");

    CtcFwdBwd c1(cNu);
    CtcFwdBwd c2(cNu2);

    CtcCompo comp(c1, c2);

    CtcFixPoint fp(comp);

    IntervalVector iV(3);
    iV[0] = Interval(NEG_INFINITY, 10);
    iV[1] = Interval(-5, 3);
    iV[2] = Interval(8, POS_INFINITY);

    fp.contract(iV);
    cout << "iV: " << iV << endl;

    //    invertEllipsoid();
    //    invertLMI();
    //    triangleBall();
    //ctcAndDraw1();
    //vibes::axisAuto();
    //ctcAndDraw2();
    //vibes::axisAuto();
    vibes::endDrawing();

    return EXIT_SUCCESS;
}

void exemple()
{
    vibes::newFigure("example");

    Matrix m1(2, 2);
    m1 = Matrix::zeros(2);
    m1[0][0] = 10;
    m1[1][1] = 1;

    Vector v1(2);
    v1[0] = 0;
    v1[0] = 0;

    double r1 = 1;

    Matrix m2(2, 2);
    m2 = Matrix::eye(2);

    Vector v2(2);
    v2[0] = 1;
    v2[1] = 0;

    double r2 = 1;

    vector<Matrix> Ps;
    Ps.push_back(m1);
    Ps.push_back(m2);

    vector<Vector> vs;
    vs.push_back(v1);
    vs.push_back(v2);

    vector<double> rs;
    rs.push_back(r1);
    rs.push_back(r2);

    CtcEllipsoids cElls(Ps, vs, rs);

    //    NumConstraint c1In("x1","x2","(x1-1)^2+x2^2>=11");
    NumConstraint c2In("x1", "x2", "x1^2+2*x1^2>=1");

    //    CtcFwdBwd cIn1(c1In);
    CtcFwdBwd cIn2(c2In);

//    CtcUnion cIn(cIn1, cIn2);

    vector<IntervalVector> in;
    vector<IntervalVector> out;
    vector<IntervalVector> perhaps;

    IntervalVector search(2, Interval(-2, 2));

//    sivia(search, cIn, cElls, in, out, perhaps);

    for (uint i = 0; i < out.size(); i++)
        vibes::drawBox(out[i], vibesParams("figure", "example", "black[blue]"));

    for (uint i = 0; i < in.size(); i++)
        vibes::drawBox(in[i], vibesParams("figure", "example", "black[red]"));

    for (uint i = 0; i < perhaps.size(); i++)
        vibes::drawBox(perhaps[i], vibesParams("figure", "example", "black[yellow]"));

}

void invertEllipsoid()
{
    vibes::newFigure("Ellipsoid inversion using inclusion test");
    vibes::setFigureProperties("Ellipsoid inversion using inclusion test", vibesParams("x", 0, "y", 500, "width", 450, "height", 450));

    vibes::newFigure("Ellipsoid inversion using Fwd-Bwd contractor");
    vibes::setFigureProperties("Ellipsoid inversion using Fwd-Bwd contractor", vibesParams("x", 460, "y", 500, "width", 450, "height", 450));

    vibes::newFigure("Ellipsoid inversion using LMI contractor");
    vibes::setFigureProperties("Ellipsoid inversion using LMI contractor", vibesParams("x", 920, "y", 500, "width", 450, "height", 450));

    Variable x(2);

    Matrix P(2, 2);
    P[0][0] = 1;
    P[1][1] = 1;
    P[0][1] = P[1][0] = 0.7;

    Vector C(2);
    C[0] = 0;
    C[1] = 0;

    double R(5);

    Function incTest(x, pow((x[0] - C[0]), 2) * P[0][0] + pow((x[1] - C[1]), 2) * P[1][1]+(x[0] - C[0])*(x[1] - C[1]) * P[0][1]+(x[0] - C[0])*(x[1] - C[1]) * P[1][0] - pow(R, 2));

    double _box[2][2] = {
        {10, 10},
        {10, 10}
    };
    IntervalVector box(2, _box);

    double _search[2][2] = {
        {-20, 20},
        {-20, 20}
    };
    IntervalVector search(2, _search);

    // SIVIA
    list<IntervalVector> l;
    l.push_back(search);

    vector<IntervalVector> out;
    vector<IntervalVector> in;
    vector<IntervalVector> perhaps;

    LargestFirst lf;

    /**
     * Test d'inclusion
     */
    for (int i = 0; i < MAX_BISSECT; i++)
    {
        IntervalVector top = l.front();
        l.pop_front();
        if (incTest.eval(top).lb() >= 0) // Exterieur
            out.push_back(top);
        else
        {
            pair<IntervalVector, IntervalVector> p = lf.bisect(top);
            l.push_back(p.first);
            l.push_back(p.second);
        }
    }

    for (uint i = 0; i < out.size(); i++)
        vibes::drawBox(out[i], vibesParams("figure", "Ellipsoid inversion using inclusion test", "black[blue]"));

    perhaps.insert(perhaps.begin(), l.begin(), l.end());

    for (uint i = 0; i < perhaps.size(); i++)
        vibes::drawBox(perhaps[i], vibesParams("figure", "Ellipsoid inversion using inclusion test", "black[yellow]"));

    Matrix P_inv(2, 2);
    real_inverse(P, P_inv);

    vibes::drawConfidenceEllipse(C[0], C[1], P_inv[0][0], P_inv[0][1], P_inv[1][1], R, vibesParams("figure", "Ellipsoid inversion using inclusion test"));


    vibes::axisAuto("Ellipsoid inversion using inclusion test");
    vibes::saveImage("Ellipsoid inversion using inclusion test", "Ellipsoid inversion using inclusion test");

    perhaps.clear();
    l.clear();
    out.clear();
    in.clear();

    /**
     * Contracteur Fwd-Bwd
     */
    NumConstraint ct(incTest, LEQ);
    NumConstraint ctIn(incTest, GEQ);
    CtcFwdBwd fwdBwd(ct);
    CtcFwdBwd fwdBwdIn(ctIn);

    sivia(search, fwdBwdIn, fwdBwd, in, out, perhaps);
    /*
        l.push_back(search);

        for (int i = 0; i < MAX_BISSECT; i++)
        {
            IntervalVector top = l.front();
            IntervalVector topBack = top;

            l.pop_front();

            try
            {
                fwdBwd.contract(top);
            }
            catch (EmptyBoxException e)
            {
                out.push_back(topBack);
            }

            if (!top.is_empty())
            {
                IntervalVector *diff;
                int n = topBack.diff(top, diff);
                for (int i = 0; i < n; i++)
                    out.push_back(diff[i]);
                pair<IntervalVector, IntervalVector> p = lf.bisect(top);
                l.push_back(p.first);
                l.push_back(p.second);
            }
        }
     */
    for (uint i = 0; i < out.size(); i++)
        vibes::drawBox(out[i], vibesParams("figure", "Ellipsoid inversion using Fwd-Bwd contractor", "black[blue]"));

    for (uint i = 0; i < in.size(); i++)
        vibes::drawBox(in[i], vibesParams("figure", "Ellipsoid inversion using Fwd-Bwd contractor", "black[red]"));


    for (uint i = 0; i < perhaps.size(); i++)
        vibes::drawBox(perhaps[i], vibesParams("figure", "Ellipsoid inversion using Fwd-Bwd contractor", "black[yellow]"));
    vibes::drawConfidenceEllipse(C[0], C[1], P_inv[0][0], P_inv[0][1], P_inv[1][1], R, vibesParams("figure", "Ellipsoid inversion using Fwd-Bwd contractor"));

    vibes::axisAuto("Ellipsoid inversion using Fwd-Bwd contractor");

    vibes::saveImage("Ellipsoid inversion using Fwd-Bwd contractor", "Ellipsoid inversion using Fwd-Bwd contractor");

    l.clear();
    out.clear();
    perhaps.clear();
    in.clear();

    /**
     * Contracteur LMI
     */

    CtcEllipsoid cEll(P_inv, C, R);

    l.push_back(search);
    sivia(search, fwdBwdIn, cEll, in, out, perhaps);

    for (uint i = 0; i < out.size(); i++)
        vibes::drawBox(out[i], vibesParams("figure", "Ellipsoid inversion using LMI contractor", "black[blue]"));

    for (uint i = 0; i < in.size(); i++)
        vibes::drawBox(in[i], vibesParams("figure", "Ellipsoid inversion using LMI contractor", "black[red]"));

    for (uint i = 0; i < perhaps.size(); i++)
        vibes::drawBox(perhaps[i], vibesParams("figure", "Ellipsoid inversion using LMI contractor", "black[yellow]"));
    vibes::drawConfidenceEllipse(C[0], C[1], P_inv[0][0], P_inv[0][1], P_inv[1][1], R, vibesParams("figure", "Ellipsoid inversion using LMI contractor"));

    vibes::axisAuto("Ellipsoid inversion using LMI contractor");
    vibes::saveImage("Ellipsoid inversion using LMI contractor", "Ellipsoid inversion using LMI contractor");
}

void invertLMI()
{
    vibes::newFigure("Naive LMI");
    vibes::setFigureProperties("Naive LMI", vibesParams("x", 0, "y", 500, "width", 450, "height", 450));

    vibes::newFigure("Non Naive LMI");
    vibes::setFigureProperties("Non Naive LMI", vibesParams("x", 0, "y", 500, "width", 450, "height", 450));

    NumConstraint n1out("x1", "x2", "(2*x1+x2*(1+sqrt(5)))/2>=0");
    NumConstraint n1in("x1", "x2", "(2*x1+x2*(1+sqrt(5)))/2<0");

    CtcFwdBwd c1out(n1out);
    CtcFwdBwd c1in(n1in);

    NumConstraint n2out("x1", "x2", "(2*x1+x2*(1-sqrt(5)))/2>=0");
    NumConstraint n2in("x1", "x2", "(2*x1+x2*(1-sqrt(5)))/2<0");

    CtcFwdBwd c2out(n2out);
    CtcFwdBwd c2in(n2in);

    CtcUnion cIn(c1in, c2in);
    CtcCompo ccOut(c1out, c2out);
    CtcFixPoint cOut(ccOut);
    IntervalVector search(2, Interval(-15, 15));

    vector<IntervalVector> out;
    vector<IntervalVector> in;
    vector<IntervalVector> perhaps;

    CtcIdentity ctcId(2);
    sivia(search, cIn, cOut, in, out, perhaps);

    for (uint i = 0; i < out.size(); i++)
        vibes::drawBox(out[i], vibesParams("figure", "Naive LMI", "black[blue]"));

    for (uint i = 0; i < in.size(); i++)
        vibes::drawBox(in[i], vibesParams("figure", "Naive LMI", "black[red]"));

    for (uint i = 0; i < perhaps.size(); i++)
        vibes::drawBox(perhaps[i], vibesParams("figure", "Naive LMI", "black[yellow]"));

    vibes::axisAuto("Naive LMI");

    /**
     * Non naif
     */
    MatrixArray mA(3, 2, 2);
    Matrix f0(2, 2);
    f0 = Matrix::zeros(2);
    Matrix f1(2, 2);
    f1 = Matrix::zeros(2); //Matrix::eye(2);
    f1[0][0] = 1;
    f1[1][1] = 1;
    Matrix f2(2, 2);
    f2 = Matrix::zeros(2);
    //f2[0][0]=0;
    //f2[1][1]=-1;
    cout << "f0: " << f0 << endl;
    cout << "f1: " << f1 << endl;
    cout << "f2: " << f2 << endl;

    f2[0][0] = 0;
    f2[1][1] = 1;
    f2[0][1] = 1;
    f2[1][0] = 1;

    mA[0] = f0;
    mA[1] = f1;
    mA[2] = f2;

    CtcBoxLMI cLMI(mA);

    search = IntervalVector(2, Interval(-15, 15));
    in.clear();
    out.clear();
    perhaps.clear();
    cout << "SIVIA" << endl;
    sivia(search, cIn, cLMI, in, out, perhaps);

    cout << "in.size(): " << in.size() << endl;
    cout << "out.size(): " << out.size() << endl;
    cout << "perhaps.size(): " << perhaps.size() << endl;

    cout << "f0: " << f0 << endl;
    cout << "f1: " << f1 << endl;
    cout << "f2: " << f2 << endl;
    for (uint i = 0; i < out.size(); i++)
        vibes::drawBox(out[i], vibesParams("figure", "Non Naive LMI", "black[blue]"));

    for (uint i = 0; i < in.size(); i++)
        vibes::drawBox(in[i], vibesParams("figure", "Non Naive LMI", "black[red]"));

    for (uint i = 0; i < perhaps.size(); i++)
        vibes::drawBox(perhaps[i], vibesParams("figure", "Non Naive LMI", "black[yellow]"));

    vibes::axisAuto("Non Naive LMI");
    vibes::saveImage("Naive LMI", "Naive LMI");
    vibes::saveImage("Non Naive LMI", "Non Naive LMI");
}

void triangleBall()
{
    double _search[2][2] = {
        {-1, 2},
        {-1, 2}
    };
    IntervalVector search(2, _search);

    vector<IntervalVector> out;
    vector<IntervalVector> in;
    vector<IntervalVector> perhaps;

    vibes::newFigure("TriangleBall FwdBwd");
    vibes::setFigureProperties("TriangleBall FwdBwd", vibesParams("x", 0, "y", 500, "width", 450, "height", 450));
    vibes::newFigure("TriangleBall LMI");
    vibes::setFigureProperties("TriangleBall LMI", vibesParams("x", 0, "y", 500, "width", 450, "height", 450));
    MatrixArray mA(3, 3, 3);

    Matrix f0(3, 3);
    f0 = Matrix::zeros(3);
    f0[0][0] = -1;
    f0[1][1] = -1;

    Matrix f1(3, 3);
    f1 = Matrix::zeros(3); //Matrix::eye(2);
    f1[0][0] = -1;
    f1[1][1] = 1;

    Matrix f2(3, 3);
    f2 = Matrix::zeros(3);
    f2[0][0] = -1;
    f2[1][1] = -1;
    f2[2][2] = 1;

    cout << "f0: " << endl << f0 << endl;
    cout << "f1: " << endl << f1 << endl;
    cout << "f2: " << endl << f2 << endl;

    mA[0] = f0;
    mA[1] = f1;
    mA[2] = f2;

    CtcBoxLMI cTOut(mA);

    NumConstraint nCTIn1("x1", "x2", "x2>=1-x1");
    NumConstraint nCTIn2("x1", "x2", "x2>=1+x1");
    NumConstraint nCTIn3("x1", "x2", "x2<=0");

    CtcFwdBwd cTIn1(nCTIn1);
    CtcFwdBwd cTIn2(nCTIn2);
    CtcFwdBwd cTIn3(nCTIn3);

    CtcUnion cTFwdBwdIn(cTIn1, cTIn2, cTIn3);


    NumConstraint nCTOut1("x1", "x2", "x2<=1-x1");
    NumConstraint nCTOut2("x1", "x2", "x2<=1+x1");
    NumConstraint nCTOut3("x1", "x2", "x2>=0");

    CtcFwdBwd cTOut1(nCTOut1);
    CtcFwdBwd cTOut2(nCTOut2);
    CtcFwdBwd cTOut3(nCTOut3);

    CtcCompo cTFwdBwdOut(cTOut1, cTOut2, cTOut3);


    MatrixArray mB(3, 3, 3);
    double det = 1 - 0.7 * 0.7;
    Matrix m0(3, 3);
    m0 = Matrix::zeros(3);
    m0[0][0] = -0.5;
    m0[1][1] = -1 / det;
    m0[2][2] = -1 / det;
    m0[1][2] = 0.7 / det;
    m0[2][1] = 0.7 / det;
    m0[1][0] = 0.5;
    m0[0][1] = 0.5;
    m0[2][0] = 0.5;
    m0[0][2] = 0.5;

    Matrix m1(3, 3);
    m1 = Matrix::zeros(3);
    m1[0][1] = 1;
    m1[1][0] = 1;

    Matrix m2(3, 3);
    m2 = Matrix::zeros(3);
    m2[0][2] = 1;
    m2[2][0] = 1;

    mB[0] = m0;
    mB[1] = m1;
    mB[2] = m2;

    CtcBoxLMI cBOut(mB);

    NumConstraint nCBIn("x1", "x2", "(x1-0.5)^2+(x2-0.5)^2+2*0.7*(x2-0.5)*(x1-0.5)>0.5");
    CtcFwdBwd cBFwdBwdIn(nCBIn);

    NumConstraint nCBOut("x1", "x2", "(x1-0.5)^2+(x2-0.5)^2+2*0.7*(x2-0.5)*(x1-0.5)<0.5");
    CtcFwdBwd cBFwdBwdOut(nCBOut);

    out.clear();
    in.clear();
    perhaps.clear();
    CtcIdentity ctcId(2);

    CtcCompo cOut(cBFwdBwdIn, cTFwdBwdOut);
    CtcCompo cOut2(cBFwdBwdOut, cTFwdBwdIn);
    CtcUnion ccOut(cOut, cOut2);

    CtcUnion cIn(cBFwdBwdOut, cTFwdBwdIn);
    CtcUnion cIn2(cBFwdBwdIn, cTFwdBwdOut);
    CtcCompo ccIn(cIn, cIn2);

    sivia(search, ccIn, ccOut/*TFwdBwdOut/*cBOut/*cOutInter*/, in, out, perhaps);

    for (uint i = 0; i < out.size(); i++)
        vibes::drawBox(out[i], vibesParams("figure", "TriangleBall FwdBwd", "black[blue]"));

    for (uint i = 0; i < in.size(); i++)
        vibes::drawBox(in[i], vibesParams("figure", "TriangleBall FwdBwd", "black[red]"));

    double area = 0;
    for (uint i = 0; i < perhaps.size(); i++)
    {
        vibes::drawBox(perhaps[i], vibesParams("figure", "TriangleBall FwdBwd", "black[yellow]"));
        area += perhaps[i].volume();
    }

    vibes::axisAuto("TriangleBall FwdBwd");
    out.clear();
    in.clear();
    perhaps.clear();

    CtcCompo cOut_bis(cBFwdBwdIn, cTOut);
    CtcCompo cOut2_bis(cBOut, cTFwdBwdIn);
    CtcUnion ccOut_bis(cOut_bis, cOut2_bis);

    CtcUnion cIn_bis(cBOut, cTFwdBwdIn);
    CtcUnion cIn2_bis(cBFwdBwdIn, cTOut);
    CtcCompo ccIn_bis(cIn_bis, cIn2_bis);

    sivia(search, ccIn_bis/*cBFwdBwdIn/*ctcId/*/, ccOut_bis/*TFwdBwdOut/*cOutInter*/, in, out, perhaps);

    for (uint i = 0; i < out.size(); i++)
        vibes::drawBox(out[i], vibesParams("figure", "TriangleBall LMI", "black[blue]"));

    for (uint i = 0; i < in.size(); i++)
        vibes::drawBox(in[i], vibesParams("figure", "TriangleBall LMI", "black[red]"));

    cout << "area FwdBwd: " << area << endl;

    area = 0;

    for (uint i = 0; i < perhaps.size(); i++)
    {
        area += perhaps[i].volume();
        vibes::drawBox(perhaps[i], vibesParams("figure", "TriangleBall LMI", "black[yellow]"));
    }

    vibes::axisAuto("TriangleBall LMI");
    cout << "area LMI: " << area << endl;
    vibes::saveImage("TriangleBall FwdBwd", "TriangleBall FwdBwd");
    vibes::saveImage("TriangleBall LMI", "TriangleBall LMI");
}

void ctcAndDraw1()
{
    vibes::newFigure("Box vs. Ellipse");
    Matrix P = Matrix::eye(2);
    P[0][0] = 0.9;
    P[0][1] = 0.8;
    P[1][0] = P[0][1];

    Matrix P_inv(2, 2);
    real_inverse(P, P_inv);

    Vector C = Vector::zeros(2);
    C[0] = 5;
    C[1] = 4;
    double R = 3;

    CtcEllipsoid ctcEll(P, C, R);

    Variable x(2);

    NumConstraint c(x, P_inv[0][0] * sqr(x[0] - C[0]) + P_inv[1][1] * sqr(x[1] - C[1]) + 2 * P_inv[0][1] * (x[0] - C[0]) * (x[1] - C[1]) <= sqr(R));
    CtcFwdBwd ctcFwdBwd(c);

    double _box[2][2] = {
        {-8, 17},
        {-2.5, 5}
    };
    IntervalVector box(2, _box);

    vibes::drawConfidenceEllipse(C[0], C[1], P[0][0], P[0][1], P[1][1], R, vibesParams("figure", "Box vs. Ellipse"));
    vibes::drawBox(box, vibesParams("figure", "Box vs. Ellipse"));
    IntervalVector boxCEll = box;
    IntervalVector boxCFwd = box;
    ctcEll.contract(boxCEll);
    vibes::drawBox(boxCEll, vibesParams("figure", "Box vs. Ellipse", "red[none]"));
    ctcFwdBwd.contract(boxCFwd);
    vibes::drawBox(boxCFwd, vibesParams("figure", "Box vs. Ellipse", "green[none]"));

    cout << " ====Ctc1==== " << endl;
    cout << "Boite originale: " << endl << box << endl;
    cout << "Boite CtcEll: " << endl << boxCEll << endl;
    cout << "Boite CtcFwdBwd: " << endl << boxCFwd << endl;
}

void ctcAndDraw2()
{
    vibes::newFigure("Box vs. 2 Ellipses");
    Matrix P1 = Matrix::eye(2);
    P1[0][0] = 0.9;
    P1[0][1] = 0.8;
    P1[1][0] = P1[0][1];

    Matrix P2 = Matrix::eye(2);
    P2[0][0] = 6;
    P2[1][1] = 0.1;

    Matrix P1_inv(2, 2);
    Matrix P2_inv(2, 2);

    real_inverse(P1, P1_inv);
    real_inverse(P2, P2_inv);

    Vector C1 = Vector::zeros(2);
    C1[0] = 5;
    C1[1] = 0;

    Vector C2 = Vector::zeros(2);

    double R1 = 3;
    double R2 = 4.5;

    vector<Matrix> Ps;
    vector<Vector> Cs;
    vector<double> Rs;

    Ps.push_back(P1);
    Ps.push_back(P2);

    Cs.push_back(C1);
    Cs.push_back(C2);

    Rs.push_back(R1);
    Rs.push_back(R2);

    CtcEllipsoids ctcEll(Ps, Cs, Rs);

    Variable x(2);

    NumConstraint c1(x, P1_inv[0][0] * sqr(x[0] - C1[0]) + P1_inv[1][1] * sqr(x[1] - C1[1]) + 2 * P1_inv[0][1] * (x[0] - C1[0]) * (x[1] - C1[1]) <= sqr(R1));
    NumConstraint c2(x, P2_inv[0][0] * sqr(x[0] - C2[0]) + P2_inv[1][1] * sqr(x[1] - C2[1]) + 2 * P2_inv[0][1] * (x[0] - C2[0]) * (x[1] - C2[1]) <= sqr(R2));

    CtcFwdBwd ctcFwdBwd1(c1);
    CtcFwdBwd ctcFwdBwd2(c2);

    CtcCompo ctcFwdBwd(ctcFwdBwd1, ctcFwdBwd2);

    double _box[2][2] = {
        {-8, 17},
        {-2.5, 5}
    };

    IntervalVector box(2, _box);

    vibes::drawConfidenceEllipse(C1[0], C1[1], P1[0][0], P1[0][1], P1[1][1], R1, vibesParams("figure", "Box vs. 2 Ellipses"));
    vibes::drawConfidenceEllipse(C2[0], C2[1], P2[0][0], P2[0][1], P2[1][1], R2, vibesParams("figure", "Box vs. 2 Ellipses"));
    vibes::drawBox(box, vibesParams("figure", "Box vs. 2 Ellipses"));

    IntervalVector boxCEll = box;
    IntervalVector boxCFwd = box;

    ctcEll.contract(boxCEll);
    vibes::drawBox(boxCEll, vibesParams("figure", "Box vs. 2 Ellipses", "red[none]"));
    ctcFwdBwd.contract(boxCFwd);
    vibes::drawBox(boxCFwd, vibesParams("figure", "Box vs. 2 Ellipses", "green[none]"));

    cout << " ====Ctc1==== " << endl;
    cout << "Boite originale: " << endl << box << endl;
    cout << "Boite CtcEll: " << endl << boxCEll << endl;
    cout << "Boite CtcFwdBwd: " << endl << boxCFwd << endl;
}

void contractAndStore(IntervalVector &box, Ctc &contractor, vector<IntervalVector> &storage)
{
    IntervalVector boxBack = box; // Get a copy
    try
    {
        contractor.contract(box);
        if (boxBack == box) // Nothing contracted
            return;
        IntervalVector* rest;
        int n = boxBack.diff(box, rest);
        for (int i = 0; i < n; i++)
            storage.push_back(rest[i]);
        delete[] rest;
    }
    catch (ibex::EmptyBoxException e)
    {
        storage.push_back(boxBack);
    }
}

void sivia(IntervalVector search, Ctc &cIn, Ctc &cOut, vector<IntervalVector> &in, vector<IntervalVector> &out, vector<IntervalVector> &perhaps)
{
    in.clear();
    out.clear();
    perhaps.clear();

    list<IntervalVector> l;
    LargestFirst lf;

    l.push_back(search);
    int nb_bis = 0;
    while (!l.empty())//for (int i = 0; i < MAX_BISSECT; i++)
    {
        IntervalVector top = l.front();

        l.pop_front();

        contractAndStore(top, cOut, out);
        contractAndStore(top, cIn, in);

        if (!top.is_empty())
        {
            if (nb_bis < MAX_BISSECT)
            {
                nb_bis++;
                pair<IntervalVector, IntervalVector> p = lf.bisect(top);
                l.push_back(p.first);
                l.push_back(p.second);
            }
            else
            {
                perhaps.push_back(top);
            }
        }
    }

    //   perhaps.insert(perhaps.begin(), l.begin(), l.end());
}
