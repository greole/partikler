#include "SPHDatastructures.hpp"
#include <algorithm>    // std::min

std::vector<float> add(
    const float &a,
    const std::vector<float> &b) {

    // todo test performance against std::transform with lambdas
    std::vector<float> ret(b.size());
    for (size_t ctr = 0; ctr < b.size(); ctr++) {
       ret[ctr] = a + b[ctr];
    }

    return ret;
}


std::vector<Vector> add(
    const std::vector<Vector> &a,
    const std::vector<Vector> &b) {

    assert(a.size() == b.size());

    // TODO test performance against std::transform with lambdas
    std::vector<Vector> ret(b.size());
    for (size_t ctr = 0; ctr < b.size(); ctr++) {
      for (size_t j = 0; j < 3; j++) {
       ret[ctr][j] = a[ctr][j] + b[ctr][j];
    }
      }

    return ret;
}

// std::vector<Point> add(
//     const std::vector<Point> &a,
//     const std::vector<Vector> &b) {

//     // TODO test performance against std::transform with lambdas
//     std::vector<Point> ret(b.size());
//     for (size_t ctr = 0; ctr < b.size(); ctr++) {
//        ret[ctr] = a[ctr] + b[ctr];
//     }

//     return ret;
// }



std::vector<float> inverse(
    const std::vector<float> &a
    ) {

    // TODO test performance against std::transform with lambdas
    std::vector<float> ret(a.size());
    for (size_t ctr = 0; ctr < a.size(); ctr++) {
       ret[ctr] = 1.0/a[ctr];
    }

    return ret;
}

std::vector<float> sum(
    const float a,
    const std::vector<float> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<float> ret(b.size());
    for (size_t ctr = 0; ctr < b.size(); ctr++) {
       ret[ctr] = a + b[ctr];
    }

    return ret;
}

std::vector<float> sum(
    const std::vector<float> &a,
    const std::vector<float> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<float> ret(a.size());
    for (size_t ctr = 0; ctr < a.size(); ctr++) {
       ret[ctr] = a[ctr] + b[ctr];
    }

    return ret;
}

// std::vector<Vector> sum(
//     const std::vector<Vector> &a,
//     const std::vector<Vector> &b) {

//     // TODO test performance against std::transform with lambdas
//     std::vector<Vector> ret(a.size());
//     for (size_t ctr = 0; ctr < a.size(); ctr++) {
//        ret[ctr] = a[ctr] + b[ctr];
//     }

//     return ret;
// }

// std::vector<Vector> multiplies(
//     const std::vector<Vector> &a,
//     const std::vector<Vector> &b) {

//     // TODO test performance against std::transform with lambdas
//     std::vector<Vector> ret(a.size());
//     for (size_t ctr = 0; ctr < a.size(); ctr++) {
//        ret[ctr] = {
//         a[ctr].x() * b[ctr].x(),
//         a[ctr].y() * b[ctr].y(),
//         a[ctr].z() * b[ctr].z()};

//     }

//     return ret;
// }

// std::vector<float> multiply(
//     const std::vector<Vector> &a,
//     const std::vector<Vector> &b) {

//     // TODO test performance against std::transform with lambdas
//     std::vector<float> ret(a.size());
//     for (size_t ctr = 0; ctr < a.size(); ctr++) {
//        ret[ctr] = a[ctr] * b[ctr];
//     }

//     return ret;
// }

// std::vector<Vector> multiply(
//     const std::vector<float> &a,
//     const std::vector<Vector> &b) {

//     // TODO test performance against std::transform with lambdas
//     std::vector<Vector> ret(a.size());
//     for (size_t ctr = 0; ctr < a.size(); ctr++) {
//        ret[ctr] = a[ctr] * b[ctr];
//     }

//     return ret;
// }

std::vector<float> multiply(
    const std::vector<float> &a,
    const std::vector<float> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<float> ret(a.size());
    for (size_t ctr = 0; ctr < a.size(); ctr++) {
       ret[ctr] = a[ctr] * b[ctr];
    }

    return ret;
}

std::vector<float> multiply(
    const float &a,
    const std::vector<float> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<float> ret(b.size());
    for (size_t ctr = 0; ctr < b.size(); ctr++) {
       ret[ctr] = a * b[ctr];
    }

    return ret;
}


// std::vector<Vector> multiply(
//     const float a,
//     const std::vector<Vector> &b) {

//     // TODO test performance against std::transform with lambdas
//     std::vector<Vector> ret(b.size());
//     for (size_t ctr = 0; ctr < b.size(); ctr++) {
//        ret[ctr] = a * b[ctr];
//     }
 
//     return ret;
// }

std::vector<float> sqr(
    const std::vector<float> &a
    ) {
    return multiply(a, a);
}

std::vector<float> divides(
    const std::vector<float> &a,
    const float b) {

    // TODO test performance against std::transform with lambdas
    std::vector<float> ret(a.size());
    for (size_t ctr = 0; ctr < a.size(); ctr++) {
       ret[ctr] = a[ctr] / b;
    }

    return ret;
}

std::vector<float> divides(
    const std::vector<float> &a,
    const std::vector<float> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<float> ret(a.size());
    for (size_t ctr = 0; ctr < a.size(); ctr++) {
       ret[ctr] = a[ctr] / b[ctr];
    }

    return ret;
}


std::vector<float> power(
    const std::vector<float> &a,
    const float b) {

    // TODO test performance against std::transform with lambdas
    std::vector<float> ret(a.size());
    for (size_t ctr = 0; ctr < a.size(); ctr++) {
       ret[ctr] = std::pow(a[ctr], b);
    }

    return ret;
}

void compute_kernel(
    Logger logger,
    const float h,
    const SPHPointField &particles,
    const SortedNeighbours &particle_neighbours,
    Kernel &kernel) {
    logger.set_scope("Kernel");
    logger.info_begin() << "Computing Kernel";
    const size_t size {particle_neighbours.ownId.size()};
    const float invh = 1 / h;
    const float W_fak2 = 21. / (256. * M_PI * h * h * h);
    const float dW_fak2 = 21. / (256. * M_PI * h * h * h * h);

    for (size_t pid = 0; pid < size; pid++) {

        // std::cout << "[DEBUG] pid" << std::endl;
        const Point &opos = particles[particle_neighbours.ownId[pid]];
        const Point &npos = particles[particle_neighbours.neighId[pid]];
        const CGALVector lenV = npos - opos;
        const float len = std::sqrt(squared_distance(opos, npos));

        const float q {len * invh};

        if (q > 2.) {
            std::cout << "[DEBUG] outside kernel radius" << std::endl;

            kernel.W[pid] = 0.0;
            kernel.dWdx[pid] = Vector {0.0, 0.0, 0.0};
            continue;
        }

        const float q3 = (q - 2.);
        const float qfac2 = q3 * q3;
        const float qfac4 = qfac2 * qfac2;

        float q2 = 2. * q;
        q2 += 1.;

        kernel.W[pid] = qfac4 * q2 * W_fak2;

        const float prefact = 10. * qfac2 * q * dW_fak2;
        if(! len == 0) {
          for(int j=0;j<3;j++) {
            kernel.dWdx[pid][j] = lenV[j]/len * prefact;
          }
        } else {
          for(int j=0;j<3;j++) {
            kernel.dWdx[pid][j] = 0.0;
          }
          std::cout << "computes.cpp:275 len == 0"  << std::endl;
        }
    }

    logger.info_end();
}

void compute_pressure_gradient(
        Logger logger,
        const SortedNeighbours &particle_neighbours,
        const Kernel             &kernel,
        const SPHFloatField &rho,
        const SPHFloatField &p,
        SPHVectorField      &dp
        ) {

    // TODO add particle mass
    // std::vector<Vector> tmp(rho.size(), {0, 0, 0});

    // SPHFloatField p_ab = p.add_ab(particle_neighbours);
    // SPHFloatField rho_ab = rho.mult_ab(particle_neighbours);

    SPHFloatField prho = p/rho;

    SPHFloatField tmp_ab = prho.add_ab(particle_neighbours);

    size_t size = particle_neighbours.ownId.size();

    // weighted sum
    for (size_t i = 0; i < size; i++) {
        size_t ownId = particle_neighbours.ownId[i];
        size_t neighId = particle_neighbours.neighId[i];

        const float val = tmp_ab[i];
        for (int j = 0; j < 3; j++) {
          dp[ownId][j]   += val * kernel.dWdx[i][j];
          dp[neighId][j] -= val * kernel.dWdx[i][j];
        }
    }

    std::cout << "dp" << dp << std::endl;

    // particle_neighbours.sum(

    //             multiply(multiply(rrho_b, p_ab), kernel.dWdx)
    //             , tmp);



    // tmp.sum()

    // dp = multiply(rrho, tmp);
}

void compute_pressure(
        Logger logger,
        const SortedNeighbours &particle_neighbours,
        const SPHFloatField &rho,
        SPHFloatField &pressure
        ) {

    const float c = 10.0;
    const float rho_0 = 1.0;
    const float gamma = 1.4;
    const float prefac = c*c*rho_0/gamma;
    const float p_0 = 1000;
    const float n1 = -1.0;

    const SPHScalarField tmp0 = (rho/rho_0).pow(gamma);
    const SPHScalarField tmp1 = ((tmp0 - 1.0)*prefac)+p_0;

    pressure.set(tmp1);
}

void compute_dnu(
        Logger logger,
        const SortedNeighbours& particle_neighbours,
        const Kernel & kernel,
        const std::vector<Vector> &velocities,
        const std::vector<float> &rho,
        std::vector<Vector> &dnu,
        const float h
        ) {
    // return cp.div(
    //         cp.filter(cp.scalar_product(dist, vel_ab), h, 0.0000001),
    //         cp.scalar_product(dist, dist))

    // logger.info_begin() << "Computing viscous term gradient";
    // std::vector<Vector> velocity_sub_ab = particle_neighbours.subtract(velocities);
    // // std::vector<Vector> velocity_sum_ab = particle_neighbours.add(velocities);
    // std::vector<float>  rho_sum_ab = particle_neighbours.add(rho);

    // std::vector<float>  v_scalar_r_ab = multiply(
    //         velocity_sub_ab,
    //         multiply(-1.0, particle_neighbours.distances));

    // std::vector<float> tmp0 = divides(v_scalar_r_ab,
    //        sum(0.001*h, particle_neighbours.squared_length));

    // std::vector<Vector> tmp1 = multiply(tmp0, kernel.dWdx);
    // std::vector<Vector> tmp2 = multiply(1.0, tmp1);
    // std::vector<Vector> tmp3 = multiply(inverse(rho_sum_ab), tmp2);

    // particle_neighbours.sum(multiply(10.0, tmp3), dnu);
    // logger.info_end();
}

void compute_density(
        Logger logger,
        const SortedNeighbours& particle_neighbours,
        const Kernel & kernel,
        SPHFloatField& rho
        ) {
    logger.info_begin() << "Computing density";
    const float zero = 0.0;
    rho.set_uniform(zero);
    rho.weighted_sum(particle_neighbours, kernel.W);

    // TODO DIRTY FIX to handle holes
    rho.lower_limit(0.001);

    logger.info_end();
}

void compute_du(
    Logger logger,
    const SortedNeighbours &particle_neighbours,
    const Kernel &kernel,
    const SPHVectorField &u,
    const SPHFloatField &rho,
    const SPHFloatField &p,
    // std::vector<Vector>& dnu,
    const SPHVectorField &dp,
    SPHVectorField &du) {
    logger.info_begin() << "Computing du/dt";
    // particle_neighbours.subtract(divides(pressure, sqr(densities)));
    // SPHFloatField pdd = (p / rho.pow(2.0)).sub_ab(particle_neighbours);
    // // std::vector<Vector> temp = sum(pdd, nu);

    // weighted sum
    // size_t size = particle_neighbours.ownId.size();
    // for (size_t i = 0; i < size; i++) {
    //     size_t ownId = particle_neighbours.ownId[i];
    //     size_t neighId = particle_neighbours.neighId[i];

    //     for (int j = 0; j < 3; j++) {
    //         const float val = dp.get_field()[i][j];
    //         du[ownId][j] -= val; // * kernel.dWdx[i][j];
    //         // du[neighId][j] += val;// * kernel.dWdx[i][j];
    //     }
    // }
  size_t size = du.size();
  std::cout << "computes.cpp du size" << size << " dp size " << dp.size() << std::endl;
  for (size_t i=0;i<size;i++) {
    for (size_t j=0;j<3;j++) {
      du[i][j] = -dp[i][j];
    }
  }
  logger.info_end();

}

void compute_u(
    Logger logger,
    const SPHVectorField &du,
    const float dt,
    SPHVectorField &u) {
  std::cout << "[DEBUG] u" << u << " du " << du << std::endl;
    u.set(u + (du * dt));
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

CGALVector facet_normal(Facet f)  {
  HalfedgeConstHandle h  = f.halfedge();
  Point              p1 = h->vertex()->point();
  Point              p2 = h->next()->vertex()->point();
  Point              p3 = h->next()->next()->vertex()->point();

  std::cout
    << "[DEBUG SLIDE]"
    << " nA " << p1
    << " nB " << p2
    << " nC " << p3
    << std::endl;
  CGALVector             n  = CGAL::cross_product(p2-p1, p3-p1);
  return n / std::sqrt(n*n);
}

struct SurfaceSlide {

  float rem; // Fraction how much of the original vector distance remains
  Point P;
  CGALVector vec; // Slide vector on current facet
  CGALVector dx; // Slide vector on current facet
  Facet_handle new_facet; // new facet of particle
  CGALVector new_facet_N;
};

SurfaceSlide slide_surface(
    // size_t facet_id,
    Point O,
    Point P,
    CGALVector PN, // original slide vector
    CGALVector D,  // face normal
    float maxDx,
    Facet_handle f) {

    // Projection of the motion vector on to the surface
    CGALVector S = PN - (D * PN) * D;

    std::cout
      << "[DEBUG SLIDE]"
      << " S " << S
      << " P " << P
      << " PN " << PN
      << std::endl;

    CGALVector OF {0, 0, 0};

    // CGALVector OP = P - O;

    // auto f = fs[facet_id];
    float frac = 1;

    // int edge_ctr = 0;

    // check if particle is already on a edge and
    // the direction of the slide vector S needs
    // to be handled
    // bool slides_on_edge = false;
    float lambda = 1.0;

    // If particle slides to next facet the remainder of dx on
    // next facet is needed
    float remainder = 0.0; 
    // for (auto edge : edges) {
    //   Point A = edge.first;
    //   Point B = edge.second;
    //   CGALVector PA = A - P;
    //   CGALVector AB = B - A;
    //   float a = AB.squared_length();
    //   float lambda1 = -(PA*AB)/a;
    //   Point X = A + lambda1*AB;


    //   // distance from intersection point X to N'
    //   CGALVector PX = X-P;
    //     //  particle is close to edge
    //     if (PX.squared_length() < 0.001) {
    //       std::cout << " SNAP to EDGE ";
    //       // correct direction
    //       float dir = AB*S;
    //       S = ((dir)/AB.squared_length())*AB;

    //       // dont slide over beginning or end;
    //       // slide direction is AB, hence towards B
    //       float dist = 0.0;
    //       if (dir >= 0) {
    //         dist = (B - P).squared_length();
    //         std::cout << " B is closest ";
    //       } else {
    //         dist = (A - P).squared_length();
    //         std::cout << " A is closest ";
    //       }
    //       std::cout << " dist " << dist;

    //       // TODO why cant maxDist be used directly?
    //       float maxDist = 0.1;

    //       if (dist < maxDist) {
    //         std::cout << " freeze to vertex ";
    //         lambda = 0.0;
    //       } else {
    //         const float oldLambda = lambda;
    //         const float newLambda = dist*0.999/S.squared_length();
    //         lambda = min(oldLambda, newLambda);
    //         std::cout << " newLambda " << newLambda << " lambda " << lambda;
    //       }

    //       S = lambda*S;

    //       slides_on_edge = true;

    //     }

    //     std::cout << " S " << S << " P+S " << P+S << std::endl;
    //     edge_ctr++;
    // }

    // if (slides_on_edge) return S;

    // check if particle crosses any edges

    auto h = Facet(*f).facet_begin();

    SurfaceSlide surface_slide {0.0, P, S, PN, f, D};

    do {
      Point A = h->vertex()->point();
      Point B = h->next()->vertex()->point();

      CGALVector PA = A - P;
      CGALVector PB = B - P;
      CGALVector AP = P - A;
      CGALVector AB = B - A;
      float a = AB.squared_length();
      float b = 2.0*AB*PA;
      float bb = b*b;
      float c = (PA).squared_length();
      float fac = 4.0*a*c;

      std::cout
        << "[DEBUG SLIDE]"
        << " A " << A
        << " B " << B
        << std::endl;

      // particle is on a vertex freeze
      // if ((c < 1.0e-02) || ( PB.squared_length()  < 1.0e-02)) {
      //   std::cout << "freeze particle on vertex" << std::endl;
      //   return S*0.0;
      // }

      // if (bb<=fac) {
      //   // std::cout << "bb " << bb << " fac " << fac << std::endl; 

      //   // std::cout << "freeze particle outside" << std::endl;
      //   // return S*0;
      //   continue;

      // };

      Point N = P + S;

      // Shortest, normal vector to facet edge, PX
      // to identify if particle move towards given
      // edge

      // intersection point X
      // from projection of P normal to AB
      // PX*AB = 0
      float lambda1 = -(PA*AB)/a;
      // Point XN is the normal intersection Point
      Point XN = A + lambda1*AB;

      // distance from intersection point X to N'
      CGALVector PXN = XN-P;

      // intersection fraction should be 0 < lambda < 1
      // float lambda1 = (-b + std::sqrt(b * b - 4 * a * c)) / (2.0 * a);
      // float len_sqr_PX = PX.squared_length();

      // // TODO this is not true for skewed triangles
      // if (lambda1 >= 1.0)
      //   {
      //     std::cout << "\tout of bounds Pos" << lambda1;
      //     // consider only motion along egde
      //     CGALVector SPr = S-((AB*S)*1.0001/AB.squared_length())*AB;
      //     float lambda2 = min(1.0, (P - B).squared_length()/SPr.squared_length());

      //     std::cout << " AB*S"<< AB*S
      //               << "SPr " << SPr
      //               << " N " <<  P + SPr
      //               << " (P-B) " << (P-B)
      //               << " |(P-B)| " << (P-B).squared_length()
      //               << " |(P-A)| " << (P-A).squared_length()
      //               << " |(A-B)| " << (A-B).squared_length()
      //               << " |SPr| " << SPr.squared_length()
      //               << " lambda2 " <<  lambda2
      //               << " N2 " <<  P + lambda2*SPr
      //               << std::endl;
      //     return SPr*lambda2;
      //   };
      // if (lambda1 <= 0 )
      //   {
      //     std::cout << "\tout of bounds Neg" << lambda1;
      //     // consider only motion along egde
      //     CGALVector SPr = ((AB*S)*0.9999/AB.squared_length())*AB-S;
      //     float lambda2 = min(1.0, (P - A).squared_length()/SPr.squared_length());

      //     // std::cout << "AB*S"<< AB*S
      //     //           << "SPr " << SPr
      //     //           << " N " <<  P + SPr
      //     //           << " (P-B) " << (P-B)
      //     //           << " |(P-B)| " << (P-B).squared_length()
      //     //           << " |(P-A)| " << (P-A).squared_length()
      //     //           << " |(A-B)| " << (A-B).squared_length()
      //     //           << " |SPr| " << SPr.squared_length()
      //     //           << " lambda2 " <<  lambda2
      //     //           << " N2 " <<  P + lambda2*SPr
      //     //           << std::endl;
      //     return SPr*lambda2;
      //   };

      // Check if negative slide vector and edge normal
      // have same direction, ie particle move towards
      // edge
      if( PXN*PN < 0) {
        std::cout
          << "[DEBUG SLIDE]"
          << "skipping edge, wrong direction"
          << std::endl;
        h++;
        continue;
      }

      // TODO check if S hits AB with 0<lambda<1
      auto res = intersection(Line(A, B), Line(P, P+S));
      Point X;
      CGALVector PX;
      if (! assign(X, res)) {
        std::cout
          << "[DEBUG SLIDE]"
          << " does not hit AB"
          << std::endl;
        h++;
        continue;
      }

      PX = X-P;

      float alpha = (PX*S)/PX.squared_length();

      float lambda2 = 1.0/alpha;

      std::cout
        << "[DEBUG SLIDE]"
        << " alpha " << alpha
        << " lambda2 " << lambda2
        << std::endl;
      // Something is strange if < lambda2<0
      if (lambda2<0.0) {
        h++;
        continue;
      }
      // if lambda2 > 1 at least the currently tested edge
      // does not restrict the particle but other edge might
      if (lambda2>1.0) {
        std::cout
          << "[DEBUG SLIDE]"
          << " unlimited "
          << std::endl;
        h++;
        continue;
      }

      // When facet has been changed, make sure vector points
      // inwards
      // SET other slide vector
      h--;
      h--;
      Facet_handle next_facet_handle = h->opposite()->facet();
      Facet next_facet = Facet(*next_facet_handle);
      CGALVector ON =  facet_normal(next_facet);
      h++;
      h++;
      OF = PN - (ON * PN) * ON;

      Point NPr = P + lambda2*PN;
      frac = min(frac,  lambda2);
      remainder = 1.0 - frac;

      // inwards pointing vector
      // assume that other component points always outwards
      Point nC;

      // TODO refactor
      nC = next_facet.halfedge()->opposite()->vertex()->point();
      std::cout
        << "[DEBUG SLIDE]"
        << " nC " << nC
        << std::endl;
      if (nC == A || nC == B) {
        nC = next_facet.halfedge()->next()->vertex()->point();
        std::cout
          << "[DEBUG SLIDE]"
          << " nC2 " << nC
          << std::endl;
      }
      if (nC == A || nC == B) {
        nC = next_facet.halfedge()->vertex()->point();
        std::cout
          << "[DEBUG SLIDE]"
          << " nC3 " << nC
          << std::endl;
      }

      CGALVector scaledS = S*frac;
      CGALVector remVec = scaledS - PN;

      CGALVector NAB = CGAL::cross_product(AB, ON);
      NAB = NAB/std::sqrt(NAB.squared_length());
      Point nN = P + scaledS;
      CGALVector NC = nC-N;

      if (NAB * NC < 0) {
        NAB = -NAB;
      }

      CGALVector PNin;
      if (NAB * remVec < 0) {
        remVec = remVec*NAB*NAB;
      }

      // TODO snap to new surface

      std::cout
        << "[DEBUG SLIDE]"
        << " A " << A
        << " B " << B
        << " P " << P
        << " S " << S
        << " P+S " << P+S
        << " P+S*frac " << nN
        << " lambdaOld " << lambda
        << " lambda2 " << lambda2
        << " N' " << NPr 
        << " ON " << ON
        << " NAB " << NAB
        << " nC " << nC
        << " NC " << NC
        << " OF " << OF
        << " remVec " << remVec
        << " remainder " << remainder
        << std::endl;

      surface_slide.rem = remainder;
      surface_slide.P = P+scaledS;
      surface_slide.vec = scaledS;
      surface_slide.dx = remVec;
      surface_slide.new_facet = next_facet_handle;
      surface_slide.new_facet_N = ON;
      h++;
    } while ( h != Facet(*f).facet_begin());

    std::cout
      << "[DEBUG SLIDE]"
      << " final frac " << frac
      << " final remainder " << remainder
      << " NPr"  << P + S*frac
      << " r*OF"  << remainder*OF
      << " newPoint "  << P + S*frac + remainder*OF
      << std::endl;

    return surface_slide;// + remainder*OF``;
}

SPHVectorField limit_dx(
    SPHVectorField &u,
    const float dt,
    const std::vector<Point> &opoints,
    FixedDistanceParticles &fdp,
    std::vector<Facet_handle> &facets,
    const SPHSizeTField &idx,
    const SPHPointField &pos
                        ) {
    SPHVectorField dx = u*dt;
    std::vector<Vector> &ret = dx.get_field();

    for (size_t i = 0; i < pos.size(); i++) {
        int mType = fdp.mType[i];
        Point O = pos[fdp.fixId[i]];
        Point P = pos[i];
        CGALVector PN = {dx[i][0], dx[i][1], dx[i][2]};
        CGALVector D = fdp.dir[i];
        float maxDx = fdp.maxDx[i];
        // const Facet facet = facets[fdp.facets[i]];

        CGALVector frac {0, 0, 0};

        std::cout
                  << "[DEBUG SLIDE]"
                  << "Particle Id " << idx.get_field()[i]
                  // // << " O (" << O << ")"
                  // << " P (" << P << ")"
                  // // << " maxDx " << maxDx
                  // // << " OP (" << OP << ")"
                  // // << " curDx " << std::sqrt(OP.squared_length())
                  // << " dx (" << PN << ")"
                  // << " |dx| " << std::sqrt(PN.squared_length())
                  // << " P+dx (" << P+PN << ")"
                  << std::endl;

        if (mType == 0) {
            frac = CGALVector {0, 0, 0};
        }

        if (mType == 2) {
          // std::cout << "point" << i << std::endl;

          float rem = 1.0;
          do {
             // Facet f =  
              auto surf_slide =
              slide_surface(O, P, PN, D, maxDx, fdp.facets[i]);
              std::cout
                << "[DEBUG SLIDE]"
                << " remainder " << rem
                << std::endl;
              P = surf_slide.P;
              PN = surf_slide.dx;
              rem = surf_slide.rem;
              fdp.facets[i] = surf_slide.new_facet;
              frac += surf_slide.vec;
              D = surf_slide.new_facet_N;
              fdp.dir[i] = D;
          } while (rem > 0.0);
            // std::cout << "surface " << i << " " << frac << std::endl;
        };

        if (mType == 3) {
            frac = half_sphere(O, P, PN, D, maxDx);
            // std::cout << "sphere " << i << " " << frac << std::endl;
        }
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

        std::cout
          << "[DEBUG SLIDE]"
          << "Particle Id " << idx.get_field()[i]
          << " DONE "
          << std::endl;
    }

    return SPHVectorField(ret, {"tmpx", "tmpy", "tmpz"}, "tmp");
}

void update_pos(
        Logger logger,
        SPHVectorField &u,
        const float dt,
        const std::vector<Point> opoints,
        FixedDistanceParticles &fpd,
        std::vector<Facet_handle>& facets,
        const SPHSizeTField &idx,
        SPHPointField& pos
        )
{
  logger.info_begin() << "Updating particle positions";

  std::cout << "u*dt" << u*dt << u << std::endl;
  SPHVectorField dx = limit_dx(u, dt, opoints, fpd, facets, idx, pos) ;
  // SPHVectorField dx = u*dt;
  // TODO set u to zero if dx == 0


  pos += dx;
  // TODO snap to stl surface

  std::cout << pos << std::endl;
  logger.info_end();

  // constrain position
}
