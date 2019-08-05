/*  Partikler - A general purpose framework for smoothed particle hydrodynamics
    simulations Copyright (C) 2019 Gregor Olenik

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    contact: go@hpsim.de
*/

#include "STLLimitedDx.hpp"

struct SkipEdge {
    bool skip;
    Point A;
    Point B;
};

struct SurfaceSlide {

    float rem; // Fraction how much of the original vector distance remains
    Point P;
    CGALVector vec;         // Slide vector on current facet
    CGALVector dx;          // Slide vector on current facet
    Facet_handle new_facet; // new facet of particle
    CGALVector new_facet_N;
    Vector u;
    SkipEdge ske;
};

// TODO Move to CGAL helpers
bool isInsideEdge(
    Point P,
    Point A,
    Point B,
    Facet f
    ) {
    // A point is inside a triangle if all shortest connections to edges
    // have negative scalar product with inside pointing edge normal
    EdgeNormal en = inwardPointingEdgeNormal(A, B, f);

    CGALVector X = Line(A, B).projection(P) - P;

    if (en.EN * X > 0)
    {return false;}
    else {
        return true;
    }
};

SurfaceSlide slide_surface(
    Point SP,      // Start position of particle
    CGALVector PN, // original slide vector
    Facet_handle f, // Handle to facet on which particle slides
    SkipEdge ske, // Skip the edgeHit detection for given edge after edge has
                  // been crossed
    Vector u) {

    Facet facet = Facet(*f);
    CGALVector FN = facet_normal(facet);
    auto h = facet.facet_begin();

    // First check if particle is still on its surface
    // otherwise edge hit detection fails
    // ie all vectors to triangle vertex are in same direction

    Plane p {
        h->vertex()->point(),
        h->next()->vertex()->point(),
        h->next()->next()->vertex()->point(),
    };

    float d = p.d() / std::sqrt(p.orthogonal_vector().squared_length());

    Point P;
    CGALVector PSP;

    P = p.projection(SP);
    // Use surface projected point as particle position
    // SP = P;

    // Projection of the motion vector on to the surface
    CGALVector S = surfaceProject(FN, PN);

    // std::cout << "[DEBUG SLIDE]"
    //           << " O " << O
    //           << " SP " << SP
    //           << " S " << S
    //           << std::endl;

    // frac keeps track of the minimal fraction of S a particle can travel
    // before hitting an edge
    float frac = 1.0;
    float lambda = 1.0;

    // If particle slides to next facet the remainder of dx on
    // next facet is needed
    float remainder = 0.0;

    Point upa;
    SkipEdge ske_ {false, upa, upa};

    SurfaceSlide surface_slide {0.0, P, S, PN, f, FN, u, ske_};

    int skipCtr = 0;
    // check if particle crosses any edges
    do {
        Point A = h->vertex()->point();
        Point B = h->next()->vertex()->point();

        // Skip current edge if particle has just crossed this edge and lies
        // directly on it
        if (ske.skip) {
            bool skipEdge =
                (A == ske.A && ske.B == B) || (A == ske.B && ske.A == B);
            if (skipEdge) {
                skipCtr++;
                h++;
                continue;
            }
        }

        bool inside = isInsideEdge(P, A, B, facet);

        // std::cout
        //   << "[DEBUG SLIDE]"
        //   << " A " << A
        //   << " B " << B
        //   << std::endl;

        CGALVector PB = B - P;
        CGALVector AP = P - A;

        if (!approachesEdge(A, B, P, S)) {
            // std::cout
            //   << "[DEBUG SLIDE] "
            //   << "skipping edge, wrong direction"
            //   << std::endl;
            skipCtr++;
            h++;
            continue;
        }

        HitPoint hit = approximateEdgeHit(A, B, P, S, 1e-16);

        if (!hit.hit) {
            // std::cout
            //   << "[DEBUG SLIDE]"
            //   << " does not hit AB"
            //   << std::endl;
            h++;
            continue;
        }

        CGALVector PX;
        PX = hit.X - P;

        float lambda = 1.0 / ((PX * S) / PX.squared_length());

        // std::cout
        //   << "[DEBUG SLIDE]"
        //   << " hit.X " << hit.X
        //   << " lambda " << lambda
        //   << std::endl;

        // if lambda > 1 at least the currently tested edge
        // does not restrict the particle but other edge might
        if (lambda > 1.0 || lambda < 0.0 || !inside) {
            // std::cout
            //   << "[DEBUG SLIDE]"
            //   << " unlimited"
            //   << " lambda " << lambda
            //   << std::endl;
            h++;
            continue;
        }

        if (lambda >= frac) {
            // std::cout
            //   << "[DEBUG SLIDE] skipping facet change"
            //   << std::endl;
            h++;
            continue;
        } else {
            frac = lambda;
        }

        remainder = 1.0 - frac;

        // When facet has been changed, make sure vector points
        // inwards
        // SET other slide vector
        h--;
        h--;
        Facet_handle next_facet_handle = h->opposite()->facet();
        Facet next_facet = Facet(*next_facet_handle);
        h++;
        h++;

        EdgeNormal en = inwardPointingEdgeNormal(A, B, next_facet);

        CGALVector ON = en.FN;  // new edge normal
        CGALVector NAB = en.EN; // new edge inward pointing normal

        CGALVector Ss = S * frac;
        CGALVector remVec = S - Ss;

        // Rotation angles around Cartesian coordinates
        // Matrix R =  createRotationMatrix(D, ON);
        Matrix R = createRotationMatrix(next_facet, facet);
        Matrix R2 = createRotationMatrix(FN, ON);

        // std::cout
        //   << "[DEBUG SLIDE] Rotate before"
        //   << " D " << D        // Current face normal
        //   << " ON " << ON      // New face normal
        //   << " R.x " << R.x    // Rotation matrix
        //   << " R.y " << R.y
        //   << " R.z " << R.z
        //   << " R2.x " << R2.x    // Rotation matrix
        //   << " R2.y " << R2.y
        //   << " R2.z " << R2.z
        //   << " remVec " << remVec // Slide vector remainder
        //   << " u " << surface_slide.u[0]
        //   << "  " << surface_slide.u[1]
        //   << "  " << surface_slide.u[2]
        //   << std::endl;

        CGALVector remVecOrig = remVec;
        if (ON != -FN) {
            // std::cout
            //   << "[DEBUG SLIDE] Rotate "
            //   << std::endl;

            remVec = rotate(R2, remVec);
            surface_slide.u = rotate(R2, u);
        }

        // std::cout
        //   << "[DEBUG SLIDE] Rotate after"
        //   << " remVec " << remVec
        //   << " remVec2 " << rotate(R, remVecOrig)
        //   << " u " << surface_slide.u[0]
        //   << "  " << surface_slide.u[1]
        //   << "  " << surface_slide.u[2]
        //   << std::endl;

        Point N = P + Ss;

        if (NAB * remVec < 0) {
            // TODO dont project but mirror
            remVec = -remVec;
            surface_slide.u[0] = -surface_slide.u[0];
            surface_slide.u[1] = -surface_slide.u[1];
            surface_slide.u[2] = -surface_slide.u[2];

            // std::cout
            //   << "[DEBUG SLIDE] invert remVec\n"
            //   << "[DEBUG SLIDE]"
            //   << " remVec " << remVec
            //   << "\n[DEBUG SLIDE]"
            //   << " u " << surface_slide.u[0]
            //   << "  " << surface_slide.u[1]
            //   << "  " << surface_slide.u[2]
            //   << std::endl;
        }

        // std::cout
        //   << "[DEBUG SLIDE] Change facet"
        //   << " N " << N
        //   << "\n[DEBUG SLIDE]"
        //   << " NAB " << NAB
        //   << "\n[DEBUG SLIDE]"
        //   << " u " << surface_slide.u[0]
        //   << "  " << surface_slide.u[1]
        //   << "  " << surface_slide.u[2]
        //   << " remVec " << remVec
        //   << " remainder " << remainder
        //   << std::endl;

        surface_slide.rem = remainder;
        surface_slide.P = P + Ss;          // + 1e-09*NAB;
        surface_slide.vec = Ss + (P - SP); // + 1e-09*NAB;
        surface_slide.dx = remVec;
        surface_slide.new_facet = next_facet_handle;
        surface_slide.new_facet_N = ON;
        surface_slide.ske = SkipEdge {true, A, B};
        h++;
    } while (h != Facet(*f).facet_begin());

    // std::cout
    //   << "[DEBUG SLIDE]"
    //   << " final frac " << frac
    //   << " final remainder " << remainder
    //   << " N "  << P + S*frac
    //   << std::endl;

    return surface_slide;
}
SPHVectorField STL_limited_dx(
    SPHVectorField &u,
    float dt,
    SPHField<Facet_handle> &facets,
    const SPHIntField &type,
    const SPHSizeTField &idx,
    const SPHPointField &pos
    ) {

    SPHVectorField dx = u*dt;
    // TODO relax dt by CFL criterion

    std::vector<Vector> &ret = dx.get_field();

#pragma omp for nowait
    for (size_t i = 0; i < pos.size(); i++) {
        int mType = type[i];
        // Point O = pos[fdp.fixId[i]];
        Point P = pos[i];
        CGALVector PN = {dx[i][0], dx[i][1], dx[i][2]};
        // CGALVector D = fdp.dir[i];
        // float maxDx = fdp.maxDx[i];

        CGALVector frac {0, 0, 0};

        // std::cout
        //           << "[DEBUG SLIDE]"
        //           << " Particle Id " << idx.get_field()[i]
        //           << std::endl;

        if (mType == 0) {
            frac = CGALVector {0, 0, 0};
        }

        if (mType == 2) {

          float rem = 1.0;
          int ctr=0;
          Point upa;
          SkipEdge ske {false, upa, upa};
          do {
              auto surf_slide =
                slide_surface(P, PN, facets[i], ske, u[i]);
              P = surf_slide.P;
              PN = surf_slide.dx;
              rem = surf_slide.rem;
              // std::cout
              //   << "[DEBUG SLIDE]"
              //   << " remainder " << rem
              //   << "\n[DEBUG SLIDE]"
              //   << " u " << u[i][0]
              //   << "  " <<  u[i][1]
              //   << "  " <<  u[i][2]
              //   << std::endl;
              facets[i] = surf_slide.new_facet;
              frac += surf_slide.vec;
              // D = surf_slide.new_facet_N;
              // fdp.dir[i] = D;
              u[i] = surf_slide.u;
              ske = surf_slide.ske;
              ctr++;
              // Hard break after 10 iterations
              // TODO use a distance/fraction based check
              if (ctr>100) break;
          } while (rem > 0.0);
        };

        // if (mType == 3) {
        //     frac = half_sphere(O, P, PN, D, maxDx);
        // }

        // Set correct velocity for fixed particles
        // TODO modify velocity for face modified particles
        if (frac[0]==0) {
          u[i][0] = 0.0;
        }
        if (frac[1]==0) {
          u[i][1] = 0.0;
        }
        if (frac[2]==0) {
          u[i][2] = 0.0;
        }

        ret[i][0] = frac.x();
        ret[i][1] = frac.y();
        ret[i][2] = frac.z();

        // std::cout
        //   << "[DEBUG SLIDE]"
        //   << " Particle Id " << idx.get_field()[i]
        //   << " DONE "
        //   << std::endl;
    }

    return SPHVectorField (ret, {"tmpx", "tmpy", "tmpz"}, "tmp", true);
}
