#include "SPHDatastructures.hpp"
#include <algorithm>    // std::min
#include <math.h>    // std::min
#include <time.h>    // std::min
#include <stdlib.h>    // std::min

void compute_kernel(
    Logger logger,
    const float h,
    const SPHPointField &particles,
    const SortedNeighbours &particle_neighbours,
    const std::vector<Facet_handle> & facets,
    Kernel &kernel) {
    logger.set_scope("Kernel");
    logger.info_begin() << "Computing Kernel";

    const size_t size {particle_neighbours.ids.size()};
    const float invh = 1.0 / h;
    // // 3d
    // const float W_fak2 = 21. / (256. * M_PI * h * h * h);
    // const float dW_fak2 = 21. / (256. * M_PI * h * h * h * h);
    // 2d
    const float W_fak2 = 7. / (64. * M_PI * h * h);
    const float dW_fak2 = 7. / (64. * M_PI * h * h * h);

    // Resize kernel
    kernel.W = std::vector<float> (size);
    kernel.dWdxo = std::vector<Vector> (size);
    kernel.dWdxn = std::vector<Vector> (size);

    for (size_t pid = 0; pid < size; pid++) {

        size_t oid = particle_neighbours.ids[pid].ownId;
        size_t nid = particle_neighbours.ids[pid].neighId;

        const Point &opos = particles[oid];
        const Point &npos = particles[nid];

        Facet F1 = Facet(*facets[oid]);
        Facet F2 = Facet(*facets[nid]);

        CGALVector lenVo = particle_neighbours.dist[pid].on;
        CGALVector lenVn = particle_neighbours.dist[pid].no;
        float len = particle_neighbours.dist[pid].len;

        const float q {len * invh};

        if (q > 2.) {
            std::cout << "computes.cpp:49 outside kernel radius"  << std::endl;
            kernel.W[pid] = 0.0;
            kernel.dWdxo[pid] = Vector {0.0, 0.0, 0.0};
            kernel.dWdxn[pid] = Vector {0.0, 0.0, 0.0};
            continue;
        }

        const float q3 = (q - 2.);
        const float qfac2 = q3 * q3;
        const float qfac4 = qfac2 * qfac2;

        float q2 = 2. * q;
        q2 += 1.;

        kernel.W[pid] = qfac4 * q2 * W_fak2;

        const float prefact = 10. * qfac2 * q * dW_fak2;
        if( len != 0.0) {
          for(int j=0;j<3;j++) {
            kernel.dWdxo[pid][j] = lenVo[j]/len * prefact;
            kernel.dWdxn[pid][j] = lenVn[j]/len * prefact;
          }
        } else {
          for(int j=0;j<3;j++) {
            kernel.dWdxo[pid][j] = 0.0;
            kernel.dWdxn[pid][j] = 0.0;
          }
          std::cout << "computes.cpp:275 len == 0"  << std::endl;
        }
    }

    logger.info_end();
}

void compute_forces(
    Logger logger,
    const SortedNeighbours &particle_neighbours,
    const SPHPointField &particles,
    const std::vector<Facet_handle> & facets,
    const float dx,
    SPHVectorField &f)
{
  // TODO runtime.submodels("force").calculate()
  // Lennard Jones Potential
  logger.info_begin() << "Computing forces";
  f.set(zeroVec);

  const size_t size {particle_neighbours.ids.size()};
  const float sigma = dx;
  const float epsilon = 1e-06;

  for (size_t pid = 0; pid < size; pid++) {

    size_t oid = particle_neighbours.ids[pid].ownId;
    size_t nid = particle_neighbours.ids[pid].neighId;
    const Point &opos = particles[oid];
    const Point &npos = particles[nid];

    Facet F1 = Facet(*facets[oid]);
    Facet F2 = Facet(*facets[nid]);

    CGALVector lenVo = particle_neighbours.dist[pid].on;
    CGALVector lenVn = particle_neighbours.dist[pid].no;

    double len = (double) particle_neighbours.dist[pid].len;

    double rat = min(4.0, max(0.01, sigma/len));
    double fi = 4.*epsilon*(pow(rat, 12.) - pow(rat, 6.) );

      for(int j=0;j<3;j++) {
        if (lenVo[j] != 0.0) {
          f[oid][j] -= (float) fi/lenVo[j];
        }
        if (lenVn[j] != 0.0) {
          f[nid][j] -= (float) fi/lenVn[j];
        }
    }
  }
  logger.info_end();
};

void add_random_noise(Logger logger, SPHVectorField &u, const float dx) {
    srand(time(NULL));
    FOR_ALL(u, ctr) {
      float x = ((float)(rand() % 100))/100.0 - 0.5;
      float y = ((float)(rand() % 100))/100.0 - 0.5;
      float z = ((float)(rand() % 100))/100.0 - 0.5;

        Vector uv = u[ctr];

        Vector rv {uv[0] + x * dx, uv[1] + y * dx, uv[2] + z * dx};
        u[ctr] = rv;
    }
}

void compute_pressure_gradient(
    Logger logger,
    const SortedNeighbours &particle_neighbours,
    const Kernel &kernel,
    const SPHFloatField &rho,
    const SPHFloatField &p,
    SPHVectorField &dp) {

    logger.info_begin() << "Computing pressure gradient";
    SPHFloatField prho = p/(rho*rho);

    SPHFloatField tmp_ab = prho.add_ab(particle_neighbours);

    size_t size = particle_neighbours.ids.size();

    // Reset
    dp.set(Vector {0,0,0});

    // weighted sum
    for (size_t i = 0; i < size; i++) {
        size_t ownId = particle_neighbours.ids[i].ownId;
        size_t neighId = particle_neighbours.ids[i].neighId;

        const float val = tmp_ab[i];
        for (int j = 0; j < 3; j++) {
          dp[ownId][j]   += val * kernel.dWdxo[i][j];
          dp[neighId][j] += val * kernel.dWdxn[i][j];
        }
    }
    // TODO Implement
    // dp.weighted_sum(particle_neighbours, tmp_ab);

    logger.info_end();
}

void compute_pressure(
        Logger logger,
        const SortedNeighbours &particle_neighbours,
        const SPHFloatField &rho,
        SPHFloatField &pressure
        ) {
    logger.info_begin() << "Computing pressure";

    const float c = 10.0;
    const float rho_0 = 1.0;
    const float gamma = 1.4;
    const float prefac = c*c*rho_0/gamma;
    const float p_0 = 1000;
    const float n1 = -1.0;

    const SPHScalarField tmp0 = (rho/rho_0).pow(gamma);
    const SPHScalarField tmp1 = ((tmp0 - 1.0)*prefac)+p_0;

    pressure = tmp1;
    logger.info_end();
}

void compute_dnu(
    Logger logger,
    const SortedNeighbours &pn,
    const Kernel &kernel,
    const SPHVectorField &u,
    const SPHFloatField &rho,
    const SPHFloatField &nu,
    const SPHPointField &pos,
    SPHVectorField &dnu) {
    logger.info_begin() << "Computing dnu";

    // const SPHFloatField tmp0 = nu.add_ab(pn) / rho.mult_ab(pn);

    const SPHVectorField dxp = particle_distance_vec(pos, pn);

    // const SPHFloatField tmp1 = (u.sub_ab(pn) * dxp) / dxp.norm();

    // const SPHFloatField tmp = tmp0*tmp1;

    const SPHFloatField tmp = (u.sub_ab(pn) * dxp)/dxp.norm();

    // TODO Reset 
    dnu.set(Vector {0,0,0});

    // weighted sum
    const size_t size = pn.ids.size();
    for (size_t i = 0; i < size; i++) {
        size_t ownId = pn.ids[i].ownId;
        size_t neighId = pn.ids[i].neighId;

        for (int j = 0; j < 3; j++) {
            dnu[ownId][j]   -= tmp[i]  * kernel.dWdxo[i][j];
            dnu[neighId][j] -= tmp[i]  * kernel.dWdxn[i][j];
        }
    }

    dnu *= nu;

    logger.info_end();
}

void compute_density(
    Logger logger,
    const SortedNeighbours &particle_neighbours,
    const Kernel &kernel,
    SPHFloatField &rho) {
    logger.info_begin() << "Computing density";

    rho.set_uniform(0.0);

    rho.weighted_sum(particle_neighbours, kernel.W);

    // rho.lower_limit(0.01);

    logger.info_end();
}

void compute_du(
    Logger logger,
    const SortedNeighbours &particle_neighbours,
    const Kernel &kernel,
    const SPHVectorField &u,
    const SPHVectorField &f,
    const SPHFloatField &rho,
    const SPHFloatField &p,
    const SPHVectorField &dnu,
    const SPHVectorField &dp,
    SPHVectorField &du) {

    logger.info_begin() << "Computing du/dt";
    // size_t size = du.size();
    // for (size_t i = 0; i < size; i++) {
    //     for (size_t j = 0; j < 3; j++) {
    //       // du[i][j] = -dp[i][j] + f[i][j];// + dnu[i][j];
    //       du[i][j] = f[i][j];// + dnu[i][j];
    //     }
    // }

    du = dnu - dp;
    logger.info_end();
}

void limit_dt_du(
         const SPHVectorField &du,
         const float maxDx,
         float & dt
         ) {
  const float maxCFL = 0.5; // DONT HARDCODE
  float maxDu = du.norm().get_max();
  std::cout << "maxDu " << maxDu << std::endl;
  float CFL = maxDu*dt*dt/maxDx;
  std::cout << "maxCFL " << CFL << std::endl;
  // max double timestep
  float two = 2.0;
  float change = min(two, maxCFL/CFL);
  dt = dt * change;
}


void compute_u(
    Logger logger,
    const SPHVectorField &du,
    const float dt,
    SPHVectorField &u) {
    logger.info_begin() << "Computing velocity";
  // TODO implement +=
    u += (du * dt);
    // std::cout << "[DEBUG] u" << u << " du " << du << std::endl;
    logger.info_end();
}

float limit_sphere_frac(CGALVector OP, CGALVector PN, float maxDx2) {
    float a = PN.squared_length();
    float b = 2.0 * OP * PN;
    float c = OP.squared_length() - maxDx2;

    float lambda1 = (-b + std::sqrt(b * b - 4 * a * c)) / (2.0 * a);
    float unity = 1.0;
    float zero = 0.0;
    return std::max(zero, std::min(unity, lambda1));
}

CGALVector
half_sphere(Point O, Point P, CGALVector PN, CGALVector D, float maxDx) {
    Point N = P + PN;
    CGALVector ON = N - O;

    float maxDx2 = maxDx * maxDx;

    // if new particle position is inside half sphere everything is fine
    // Positive dot product means new particle position is
    // in upper half of the sphere, as defined by extrusion dir
    bool newPosIsUpperHalfSphere = (ON * D >= 0);
    // Distance to original point is less than maxDx
    bool newPosIsInsideSphere = ON.squared_length() <= maxDx2;
    if (newPosIsUpperHalfSphere && newPosIsInsideSphere) {
        return PN;
    };

    CGALVector OP = P - O;
    float a = PN.squared_length();
    bool currentPosInsideSphere = OP.squared_length() <= maxDx2;

    // if particle is currently inside sphere but new position
    // would be outside the sphere dx is to be limited
    float frac = 0;
    if (currentPosInsideSphere && (a > 0)) {
        frac = limit_sphere_frac(OP, PN, maxDx2);

        // if dot product of extrusion dir and dx is negative
        // particle travels towards half sphere boundary thus frac
        // is potentially to be limited
        float denom = PN * D;
        bool canCrossHalfSphereBound = denom < 0;
        if (canCrossHalfSphereBound) {
            float nom = -OP * D;
            float lambda2 = nom / denom;
            float zero = 0.0;
            frac = std::max(zero, std::min(frac, lambda2));
        }
    }
    return PN * frac;
}

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

// Matrix facetCoordSys(Facet f, CGALVector S) {
//   CGALVector norm = facet_normal(f);
//   CGALVector surf = S/std::sqrt(S.squared_length());
//   CGALVector third = cross_product(surf, norm);

//   return {surf, norm, third};
// }

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

SPHVectorField limit_dx(
    SPHVectorField &u,
    float dt,
    const std::vector<Point> &opoints,
    std::vector<Facet_handle> &facets,
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


void snap_to_surface(
  std::vector<Facet_handle>& facets,
  SPHPointField& pos
                     ){

  std::vector<Vector> corrector(pos.size(), {0,0,0});

  for (size_t i = 0; i < pos.size(); i++) {
  // compute distance to facet
    // Facet facet = Facet(*facets[i]);
    // CGALVector N = facet_normal(facet);
    // float d = facet.plane().d();
    // CGALVector P = pos[i]-Point(0, 0, 0);
    // float dist = (N*P+d)/std::sqrt(N.squared_length());

    Facet facet = Facet(*facets[i]);
    auto h = facet.facet_begin();

    Plane p {
             h->vertex()->point(),
             h->next()->vertex()->point(),
             h->next()->next()->vertex()->point(),
    };

    Point P = p.projection(pos[i]);
    CGALVector dx = P - pos[i];

    // if (dist != 0.0) {
    //   std::cout
    //     << "[DEBUG SLIDE] Snap to STL"
    //     << " Particle Id " << i
    //     << " Pos " << P
    //     << " dist " << dist
    //     << " N " << N
    //     << " d*N" << d*N
    //     << " p-d*N " << pos[i] - d*N
    //     << " DONE "
    //     << std::endl;
        corrector[i] = Vector{(float)dx.x(), (float)dx.y(), (float)dx.z()};
  }

  SPHVectorField cor = SPHVectorField(corrector, {"tmpx", "tmpy", "tmpz"}, "tmp");
}

float update_pos(
        Logger logger,
        SPHVectorField &u,
        float dt,
        float dx_max,
        const std::vector<Point> opoints,
        // FixedDistanceParticles &fpd,
        const SPHIntField &type,
        std::vector<Facet_handle>& facets,
        const SPHSizeTField &idx,
        SPHPointField& pos
        )
{
  logger.info_begin() << "Updating particle positions";
  const float maxCFL = 0.5; // TODO DONT HARDCODE

  SPHPointField old_pos(pos);

  SPHVectorField dx = limit_dx(u, dt, opoints, facets, type, idx, pos) ;

  pos += dx;

  const float CFL = 0.01;
  float current_max_dx = (pos - old_pos).norm().get_max();
  float dx_ratio = CFL*dx_max/current_max_dx;
  float two = 10.0;
  float change = min(two, dx_ratio);
  dt =  dt * change;

  logger.info_end();
  return dt;
}


// SPHPointField cut(SPHPointField a, SPHPointField b, Surfaces) {
//   // Deletes all Points from a that are outside of b
  
// }
