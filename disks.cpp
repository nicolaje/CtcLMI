#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "CtcLMI.h"
#include "vibes/vibes.h"
#include "CtcEllipsoid.h"
#include "CtcEllipsoids.h"


using namespace ibex;
using namespace std;

void ctcAndDraw1();
void ctcAndDraw2();

int main()
{
    vibes::beginDrawing();
    ctcAndDraw1();
    //    ctcAndDraw2();
    vibes::endDrawing();
    return EXIT_SUCCESS;
}

// Differentes boites par rapport a un disque

void ctcAndDraw1()
{
    vibes::newFigure("Boxes vs 1 disk");
    Vector C(2);
    double R = 2;

    Variable x(2);

    NumConstraint c(x, sqr(x[0] - C[0]) + sqr(x[1] - C[1]) <= sqr(R));
    CtcFwdBwd ctcFwdBwd(c);

    vibes::drawEllipse(C[0], C[1], R, R, 0, vibesParams("figure", "Boxes vs 1 disk", "lightGray[lightGray]"));

    double _box[2][2] = {
        {-4, 4},
        {-8, 3}
    };

    IntervalVector box(2, _box);
    double _smallBox[2][2] = {
        {1, 1.8},
        {-1.8, -1}};

    IntervalVector smallBox(2, _smallBox);

    vibes::drawBox(box, vibesParams("figure", "Boxes vs 1 disk"));
    vibes::drawBox(smallBox, vibesParams("figure", "Boxes vs 1 disk"));

    IntervalVector box1 = box;
    ctcFwdBwd.contract(box1);
    vibes::drawBox(box1, vibesParams("figure", "Boxes vs 1 disk"));
    box1=smallBox;
    ctcFwdBwd.contract(box1);
    vibes::drawBox(box1, vibesParams("figure", "Boxes vs 1 disk"));
    vibes::axisAuto("Boxes vs 1 disk");
}