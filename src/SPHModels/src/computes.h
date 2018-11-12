#include "SPHDatastructures.h"

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
       ret[ctr] = a[ctr] + b[ctr];
    }

    return ret;
}

std::vector<Point> add(
    const std::vector<Point> &a,
    const std::vector<Vector> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<Point> ret(b.size());
    for (size_t ctr = 0; ctr < b.size(); ctr++) {
       ret[ctr] = a[ctr] + b[ctr];
    }

    return ret;
}

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

std::vector<Vector> sum(
    const std::vector<Vector> &a,
    const std::vector<Vector> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<Vector> ret(a.size());
    for (size_t ctr = 0; ctr < a.size(); ctr++) {
       ret[ctr] = a[ctr] + b[ctr];
    }

    return ret;
}

std::vector<Vector> multiplies(
    const std::vector<Vector> &a,
    const std::vector<Vector> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<Vector> ret(a.size());
    for (size_t ctr = 0; ctr < a.size(); ctr++) {
       ret[ctr] = {
        a[ctr].x() * b[ctr].x(),
        a[ctr].y() * b[ctr].y(),
        a[ctr].z() * b[ctr].z()};

    }

    return ret;
}

std::vector<float> multiply(
    const std::vector<Vector> &a,
    const std::vector<Vector> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<float> ret(a.size());
    for (size_t ctr = 0; ctr < a.size(); ctr++) {
       ret[ctr] = a[ctr] * b[ctr];
    }

    return ret;
}

std::vector<Vector> multiply(
    const std::vector<float> &a,
    const std::vector<Vector> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<Vector> ret(a.size());
    for (size_t ctr = 0; ctr < a.size(); ctr++) {
       ret[ctr] = a[ctr] * b[ctr];
    }

    return ret;
}

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


std::vector<Vector> multiply(
    const float a,
    const std::vector<Vector> &b) {

    // TODO test performance against std::transform with lambdas
    std::vector<Vector> ret(b.size());
    for (size_t ctr = 0; ctr < b.size(); ctr++) {
       ret[ctr] = a * b[ctr];
    }

    return ret;
}

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
        const ParticleNeighbours& particle_neighbours,
        Kernel& kernel
        )
{
    logger.set_scope("Kernel");
    logger.info_begin() << "Computing Kernel";
    const size_t size {particle_neighbours.origId.size()};
    const float invh = 1/h;
    const float W_fak2  = 21. / (256. * M_PI * h * h * h);
    const float dW_fak2 = 21. / (256. * M_PI * h * h * h * h);

    for (size_t pid = 0; pid < size; pid++) {

        const float len = std::sqrt(particle_neighbours.squared_length[pid]);

        const float q {len * invh};

        if (q > 2.) {
            std::cout << "[DEBUG] outside kernel radius" << std::endl;
            kernel.W[pid] = 0.0;
            kernel.dWdx[pid] = Vector {0.0, 0.0, 0.0};
            continue;
        }

        const float q3 = (q - 2.);
        const float qfac2 = q3*q3;
        const float qfac4 = qfac2*qfac2;

        float q2 = 2.*q;
        q2 += 1.;

        kernel.W[pid] = qfac4*q2*W_fak2;

        const float prefact = 10. * qfac2 * q * dW_fak2;
        kernel.dWdx[pid] = particle_neighbours.normalised_distances[pid]*prefact;
    }

    logger.info_end();
}

void compute_pressure_gradient(
        Logger logger,
        const ParticleNeighbours &particle_neighbours,
        const Kernel             &kernel,
        const std::vector<float> &rho,
        const std::vector<float> &p,
        std::vector<Vector>      &dp
        ) {

    // TODO add particle mass
    const std::vector<float> rrho = inverse(rho);
    const std::vector<float> rrho_b = particle_neighbours.get_b(rrho);

    std::vector<Vector> tmp(rho.size(), {0, 0, 0});

    std::vector<float> p_ab = particle_neighbours.add(p);

    particle_neighbours.sum(
                multiply(multiply(rrho_b, p_ab), kernel.dWdx)
                , tmp);

    dp = multiply(rrho, tmp);
}


void compute_pressure(
        Logger logger,
        const ParticleNeighbours& particle_neighbours,
        const std::vector<float> &rho,
        std::vector<float> &pressure
        ) {

    const float c = 10.0;
    const float rho_0 = 1.0;
    const float gamma = 1.4;
    const float prefact= c*c*rho_0/gamma;
    const float p_0 = 1000;
    const float n1 = -1.0;

    const std::vector<float> tmp0 = divides(rho, rho_0);
    const std::vector<float> tmp1 = power(tmp0, gamma);
    const std::vector<float> tmp2 = add(n1, tmp1);
    const std::vector<float> tmp3 = add(p_0, tmp2);

    pressure = multiply(prefact, tmp3);
}

void compute_dnu(
        Logger logger,
        const ParticleNeighbours& particle_neighbours,
        const Kernel & kernel,
        const std::vector<Vector> &velocities,
        const std::vector<float> &rho,
        std::vector<Vector> &dnu,
        const float h
        ) {
    // return cp.div(
    //         cp.filter(cp.scalar_product(dist, vel_ab), h, 0.0000001),
    //         cp.scalar_product(dist, dist))

    logger.info_begin() << "Computing viscous term gradient";
    std::vector<Vector> velocity_sub_ab = particle_neighbours.subtract(velocities);
    // std::vector<Vector> velocity_sum_ab = particle_neighbours.add(velocities);
    std::vector<float>  rho_sum_ab = particle_neighbours.add(rho);

    std::vector<float>  v_scalar_r_ab = multiply(
            velocity_sub_ab,
            multiply(-1.0, particle_neighbours.distances));

    std::vector<float> tmp0 = divides(v_scalar_r_ab,
           sum(0.001*h, particle_neighbours.squared_length));

    std::vector<Vector> tmp1 = multiply(tmp0, kernel.dWdx);
    std::vector<Vector> tmp2 = multiply(1.0, tmp1);
    std::vector<Vector> tmp3 = multiply(inverse(rho_sum_ab), tmp2);

    particle_neighbours.sum(multiply(10.0, tmp3), dnu);
    logger.info_end();
}

void compute_density(
        Logger logger,
        const ParticleNeighbours& particle_neighbours,
        const Kernel & kernel,
        std::vector<float>& rho
        ) {
    logger.info_begin() << "Computing density";
    std::fill(rho.begin(), rho.end(), 0);
    particle_neighbours.sum(kernel.W, rho);
    logger.info_end();
}


void compute_du(
        Logger logger,
        const ParticleNeighbours& particle_neighbours,
        const Kernel & kernel,
        std::vector<Vector>& du,
        std::vector<Vector>& velocities,
        std::vector<float>& densities,
        std::vector<Vector>& dnu,
        std::vector<Vector>& pressure_gradient,
        const float dx
        )
{
        // std::vector<float> pdd = particle_neighbours.subtract(divides(pressure, sqr(densities)));
        // std::vector<Vector> temp = sum(pdd, nu);

        // particle_neighbours.sum(multiply(temp, kernel.dWdx), du);

        du = sum(multiply(-1.0, pressure_gradient), dnu);
        // du = multiply(-1.0, pressure_gradient);
}


void compute_u(
        Logger logger,
        std::vector<int>&  type,
        std::vector<Vector>& du,
        std::vector<Vector>& velocities,
        const double dt
        )
{
  const std::vector<Vector> tmp = add(velocities, multiply(dt, du));
  velocities = tmp;

  size_t n = velocities.size();
  for (size_t i=0; i<n; i++) {
      if (type[i] == 0) { velocities[i] = Vector {0,0,0};}
  }
}
void update_pos(
        Logger logger,
        std::vector<Point>& positions,
        const std::vector<Vector>& velocities,
        const double dt
        )
{
  const std::vector<Point> tmp = add(positions, multiply(dt, velocities));
  positions = tmp;
}
