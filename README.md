# icassp-oov-recognition

This has data and code related to the paper accepted at ICASSP21 "A comparison of methods for OOV-word recognition on a new Public Dataset": 

# data

This contains for English and German:  
    - The train and test set in kaldi format (audio files not included)  
    - The lexicon  
    - For convenience the list of OOVs in the test set relative to the lexicon  
    - The lexicon for the OOV-words

English LM data: [link](http://www.mediafire.com/file/fy8841cfkwft5tu/en_lm_text.txt.gz/file)
German LM data: [link](http://www.mediafire.com/file/7egjt3mygxk6whw/de_lm_text.txt.gz/file)

# scripts

Currently contains scripts to  
    - create the train/test partition from a (kaldi formatted) data folder containing CommonVoice data, `build_cv_test_train.py`  
    - create the HCL graph which can be inserted into an existing HCLG, `compose_hcl.sh`  
    - recover words from a decoded lattice that phones arcs attached to the `<unk>` token, `recover_unk_words.sh`  

# libs

More up-to-date version of wrapper is [here](https://github.com/RuABraun/fst-util). You should use that.

---

This has code which wraps OpenFST, and functions for modifying graphs (`insert`, `replace_single`, `add_boost`).

To compile you will need to include add a symlink inside the libs/ directory to a copy of the pybind11 repository, and to use `LD_LIBRARY_PATH` needs have the OpenFST libs in its path and copy the compiled .so to the site-packages/ directory (run `python -m site` to find).

# How to add words to HCLG

As mentioned in the paper, this method requires you to use a monophone model. Additionally, your language model needs to have been trained with pocolm and the `--limit-unk-history` option.

For simplicity, the modification is done on a graph without self-loops. So you need to modify `utils/mkgraph.sh` and comment L167: `rm $dir/HCLGa.fst $dir/Ha.fst 2>/dev/null || true` because we will use `HCLGa.fst`.

Inside the graph dir where the HCLG is there is a `words.txt`. You need to assign IDs to the new words you're adding and append these to `words.txt` file (these should be larger than the existing ones obviously).

Assuming all this is ready you can use `script/compose_hcl.sh` to create the HCL from a lexicon of the OOV words you want to add. Check the script for the input arguments, `model` is the `final.mdl`, isym is phones osym words. Notice it uses `create_lfst.py` so you need to fst wrapper installed. There is one hardcoded parameter on L25, `303`, see [here](https://groups.google.com/g/kaldi-help/c/jL8VnwKGRWs/m/-Pe29-G9AgAJ) for what's about. You can set it to any number larger than the existing phone IDs.

After calling the script and creating the `HCL.fst` you use the fst wrapper to modify the `HCLGa.fst`.

```
from wrappedfst import WrappedFst
fst = WrappedFst('HCLGa.fst')
ifst = WrappedFst('HCL.fst')
unk_id =  # unk symbol
fst.replace_single(unk_id, ifst)
fst.write('HCLGa_new.fst')
```

Then add the self-loops (check `mkgraph.sh` for how to do that) and you are done. Replace an existing `HCLG.fst` with the new version and you can run decoding as you would normally.
