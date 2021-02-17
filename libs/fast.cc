// Copyright (c) 2021 Idiap Research Institute, http://www.idiap.ch/
// Written by Rudolf A. Braun <rbraun@idiap.ch>
//
// This file is part of icassp-oov-recognition
//
// icassp-oov-recognition is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 3 as
// published by the Free Software Foundation.
//
// icassp-oov-recognition is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with icassp-oov-recognition. If not, see <http://www.gnu.org/licenses/>.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include <tuple>
#include <queue>
#include <iostream>
#include "fst/fst.h"
#include "fst/script/fstscript.h"
#include "fst-wrapper.h"
#include <math.h>
#include <chrono>

namespace py = pybind11;


PYBIND11_MODULE(wrappedfst, m) {
  m.doc() = "pybind11 plugin";

  py::class_<Arc>(m, "Arc")
    .def(py::init<int, int, double, int>())
    .def_readwrite("ilabel", &Arc::ilabel)
    .def_readwrite("olabel", &Arc::olabel)
    .def_readwrite("weight", &Arc::weight)
    .def_readwrite("nextstate", &Arc::nextstate);

  py::class_<WrappedFst>(m, "WrappedFst")
    .def(py::init<>())
    .def(py::init<std::string>())
    .def("read", &WrappedFst::Read)
    .def("write", &WrappedFst::Write)
    .def("write_ark_entry", &WrappedFst::WriteArkEntry)
    .def("read_ark_entries", &WrappedFst::ReadArkEntries)
    .def("add_state", &WrappedFst::AddState)
    .def("set_start", &WrappedFst::SetStart)
    .def("set_final", &WrappedFst::SetFinal, py::arg("state"),py::arg("weight")=0.)
    .def("add_arc", &WrappedFst::AddArc)
    .def("get_start", &WrappedFst::GetStart)
    .def("get_arcs", &WrappedFst::GetArcs)
    .def("determinize", &WrappedFst::Determinize)
    .def("minimize", &WrappedFst::Minimize)
    .def("arc_sort", &WrappedFst::ArcSort)
    .def("compose", &WrappedFst::Compose)
    .def("shortest_path", &WrappedFst::ShortestPath)
    .def("final", &WrappedFst::Final)
    .def("is_final", &WrappedFst::isFinal)
    .def("states", &WrappedFst::States)
    .def("connect", &WrappedFst::Connect)
    .def("delete_arcs", &WrappedFst::DeleteArcs)
    .def("delete_states", &WrappedFst::DeleteStates)
    .def("num_states", &WrappedFst::NumStates)
    .def("num_arcs", &WrappedFst::NumArcs)
    .def("insert", &WrappedFst::Insert)
    .def("replace_single", &WrappedFst::ReplaceSingle)
    .def("add_boost", &WrappedFst::AddBoost)
    .def("normalise_weights", &WrappedFst::NormaliseWeights)
    .def("copy", &WrappedFst::Copy,  py::return_value_policy::take_ownership)
    .def(py::pickle(
      [](const WrappedFst& f) {
        int num_states = f.NumStates();
        int start_state = f.GetStart();
        py::list final_states;
        py::list arcs;
        for (int state: f.States()) {
          if (f.Final(state) == 0.) {
            final_states.append(py::int_(state));
          }
          for (Arc& arc: f.GetArcs(state)) {
            arcs.append(py::make_tuple(py::int_(state), py::int_(arc.nextstate), py::int_(arc.ilabel), py::int_(arc.olabel), py::float_(arc.weight)));
          }
        }
        return py::make_tuple(num_states, start_state, final_states, arcs);
      },
      [](py::tuple t) {
        WrappedFst f;
        int num_states = t[0].cast<int>();
        int start_state = t[1].cast<int>();
        py::list final_states = t[2];
        py::list arcs = t[3];
        for (int i = 0; i < num_states; i++) {
          int state = f.AddState();
          if (state == start_state) f.SetStart(state);
        }
        for (py::handle obj: final_states) {
          int state = obj.cast<int>();
          f.SetFinal(state);
        }
        for (py::handle obj: arcs) {
          py::tuple tpl = obj.cast<py::tuple>();
          f.AddArc(tpl[0].cast<int>(), tpl[1].cast<int>(), tpl[2].cast<int>(), tpl[3].cast<int>(), tpl[4].cast<double>());
        }
        return f;
      }
      ))
      .def("__copy__", [](const WrappedFst& wfst) {
        return WrappedFst(wfst);
      })
      .def("__deepcopy__", [](const WrappedFst& wfst) {
        return WrappedFst(wfst);
      });

  py::class_<ArcIterator>(m, "ArcIterator")
    .def(py::init<WrappedFst&, int>())
    .def("Done", &ArcIterator::Done)
    .def("Next", &ArcIterator::Next)
    .def("Value", &ArcIterator::Value)
    .def("SetValue", &ArcIterator::SetValue);

}
