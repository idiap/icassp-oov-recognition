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

#include "fst/script/fstscript.h"
#include<string>
#include<vector>


struct Arc {
  int ilabel, olabel, nextstate;
  double weight;
  Arc(int ilabel, int olabel, double weight, int nextstate): ilabel(ilabel), olabel(olabel), weight(weight), nextstate(nextstate) {}
  Arc() {}
};


class WrappedFst {
public:
  fst::script::VectorFstClass* fst_;
  WrappedFst() {
    fst_ = new fst::script::VectorFstClass("standard");
  }

  WrappedFst(std::string fst_fpath) {
    fst_ = new fst::script::VectorFstClass("standard");
    Read(fst_fpath);
  }

  WrappedFst(const WrappedFst& wfst) {
    fst_ = new fst::script::VectorFstClass("standard");
    int start_state = wfst.GetStart();
    std::vector<int> final_states;
    for (int state: wfst.States()) {
      if (wfst.Final(state) == 0.) {
        final_states.push_back(state);
      }
      this->AddState();
    }

    for (int state: wfst.States()) {
      for (Arc& arc: wfst.GetArcs(state)) {
        this->AddArc(state, arc.nextstate, arc.ilabel, arc.olabel, arc.weight);
      }
    }
    this->SetStart(start_state);
    for (int state: final_states) {
      this->SetFinal(state);
    }
  }

  int AddState();

  void SetStart(int state);

  void SetFinal(int state, double weight=0.);

  void AddArc(int start_state, int next_state, int ilabel, int olabel, double weight);

  void Read(std::string fst_fpath);

  static std::vector<std::pair<std::string, WrappedFst*>> ReadArkEntries(std::string fst_fpath) {
    std::fstream fs;
    fs.open(fst_fpath, std::fstream::binary | std::fstream::in);
    const fst::FstReadOptions opts(fst_fpath);
    std::vector<std::pair<std::string, WrappedFst*>> lst;
    while (!fs.eof()) {
      std::string key;
      fs >> key;
      if (key.empty()) break;
      int space = fs.get();
      WrappedFst* fst = new WrappedFst();
      fst->fst_ = fst::script::VectorFstClass::Read<fst::StdArc>(fs, opts);
      lst.emplace_back(key, fst);
    }
    fs.close();
    return lst;
  }

  void Write(std::string fst_fpath);

  void WriteArkEntry(std::string key, std::string fst_fpath);

  int GetStart() const;

  double Final(int state) const;

  bool isFinal(int state) const { return Final(state) != std::numeric_limits<double>::infinity(); }

  std::vector<Arc> GetArcs(int state) const;

  void Determinize();

  void Minimize();

  void ArcSort(std::string);

  void Compose(WrappedFst& other);

  void ShortestPath();

  std::vector<int> States() const;

  void Connect();

  void DeleteArcs(int state);

  void DeleteStates(std::vector<int64_t> states) { fst_->DeleteStates(states); }

  int NumStates() const;

  int NumArcs(int state) const;

  WrappedFst* Copy() const;

  void Insert(const int olabel, WrappedFst* fst);

  void ReplaceSingle(const int olabel, WrappedFst* fst);

  void AddBoost(std::vector< std::vector<int>> word_subwords, double boost, int disambig, int unk);

  void NormaliseWeights();

//  bool CheckHasEpsilonLoop(int start, int end);

  ~WrappedFst() {
    delete fst_;
  }
};

class ArcIterator {
public:
  fst::script::MutableArcIteratorClass arc_iterator;
  int num_arcs;
  int count_next_;
  ArcIterator(WrappedFst& fst, int state): arc_iterator(fst.fst_, state) {
    num_arcs = fst.fst_->NumArcs(state);
    count_next_ = 0;
  }

  ~ArcIterator() { }

  bool Done() { bool done = arc_iterator.Done(); if (done) count_next_ = 0; return done; }
  void Next() { arc_iterator.Next(); ++count_next_; }

  void NextI(int i) {
    for (int j = 0; j < i; ++j) {
      arc_iterator.Next();
    }
  }

  Arc Value() {
    fst::script::ArcClass arcc = arc_iterator.Value();
    Arc arc(arcc.ilabel, arcc.olabel, arcc.weight.GetWeight<fst::TropicalWeight>()->Value(), arcc.nextstate);
    return arc;
  }

  void SetValue(Arc arc) {
    fst::StdArc sarc(arc.ilabel, arc.olabel, fst::TropicalWeight(arc.weight), arc.nextstate);
    arc_iterator.SetValue(fst::script::ArcClass(sarc));
  }
};
