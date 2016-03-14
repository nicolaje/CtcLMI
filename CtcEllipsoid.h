#ifndef CTCELLIPSOID_H
#define	CTCELLIPSOID_H

#include <ibex/ibex.h>
#include "CtcLMI.h"

class CtcEllipsoid : public ibex::Ctc {
public:

    CtcEllipsoid(const ibex::Matrix &P, const ibex::Vector &C, const double &R) : Ctc(-1), P(P), C(C), R(R) {
    };
    void contract(ibex::IntervalVector &box);

private:
    const ibex::Matrix &P;
    const ibex::Vector &C;
    const double &R;
    void mkMatrixArray(const ibex::IntervalVector &box);
    ibex::MatrixArray *matrices;
};

#endif	/* CTCELLIPSOID_H */

