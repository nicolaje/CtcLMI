#ifndef CTCELLIPSOIDS_H
#define	CTCELLIPSOIDS_H

#include <ibex/ibex.h>
#include <vector>
#include "CtcLMI.h"

class CtcEllipsoids : public ibex::Ctc {
public:
    CtcEllipsoids(const std::vector<ibex::Matrix> &Ps, const std::vector<ibex::Vector> &Cs, const std::vector<double> &Rs): Ctc(-1), Ps(Ps), Cs(Cs), Rs(Rs){};
    void contract(ibex::IntervalVector &box);

private:
    const std::vector<ibex::Matrix> &Ps;
    const std::vector<ibex::Vector> &Cs;
    const std::vector<double> &Rs;
    void mkMatrixArray(const ibex::IntervalVector &box);
    ibex::MatrixArray *matrices;
};

#endif	/* CTCELLIPSOIDS_H */

