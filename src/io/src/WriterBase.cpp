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

#include "WriterBase.hpp"

WriterBase::WriterBase(
    const std::string &model_name, YAML::Node parameter, ObjectRegistry &objReg)
    : Model(model_name, parameter, objReg),
      write_freq_(read_or_default_coeff<int>("freq", 1.0)),
      last_write_(0), // TODO seek on disk for last write
      time_graph_(objReg.get_object<TimeGraph>("TimeGraph")) {}

bool WriterBase::write() {
    int ct = time_graph_.get_current_timestep();
    if (ct % write_freq_ == 0) return true;
    return false;
}
