(* Distance-Based Boolean Applicability Domain (DBBAD)
  reference implementation.

  Copyright (C) 2020, Francois Berenger

  Yamanishi laboratory,
  Department of Bioscience and Bioinformatics,
  Faculty of Computer Science and Systems Engineering,
  Kyushu Institute of Technology,
  680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

open Printf
open Sys

module CLI = Minicli.CLI
module A = MyArray
module FpMol = Molenc.FpMol
module Ht = BatHashtbl
module L = MyList
module Log = Dolog.Log
module Mol = Molecules

let write_new_dataset mol =
  let output = open_out "new_dataset.txt" in
  output_string output mol ^ "\n";
  close_out output

let apply_DBBAD ncores best_d dataset active_molecules =
  let actives = Molecules.from_file active_molecules in
  let actives_bst = Bstree.(create 1 Two_bands (A.of_list actives)) in
  let maybe_ok_test_mols =
    Parany.Parmap.parmap ncores (fun test_mol ->
        (Dbad_common.mol_is_inside_global_AD test_mol best_d actives_bst,
          test_mol)
      ) dataset in
  let ok_test_mols =
    L.fold (fun acc (maybe_ok, mol) ->
        if maybe_ok then begin
          write_new_dataset mol;
          mol :: acc 
        end
      ) [] maybe_ok_test_mols in
  let ok_card = L.length ok_test_mols in
  let test_card = L.length dataset in
  Log.info "passed AD: %d / %d" ok_card test_card;
  ok_test_mols

let process_file dataset_folder actives file ncores best_d nfolds dscan_fn =
  let chunk_dbbad_fn = file ^ ".dbbad" in
  let chunk = Molecules.from_file file in
  let chunk_DBBAD = apply_DBBAD ncores best_d chunk actives in
  Utls.with_out_file chunk_dbbad_fn (fun out ->
      L.iter (FpMol.to_out out) test_DBBAD
    );
  Log.info "chunk_DBBAD written to: %s" chunk_dbbad_fn;
  (* report actives proportion before/after DBBAD *)
  let card_act_before, card_dec_before =
    L.filter_counts FpMol.is_active chunk in
  let card_act_after, card_dec_test_after =
    L.filter_counts FpMol.is_active chunk_DBBAD in
  let old_rate =
    let atrain = float card_act_before in
    let dtrain = float card_dec_before in
    atrain /. (atrain +. dtrain) in
  let new_rate =
    let atest = float card_act_after in
    let dtest = float card_dec_test_after in
    atest /. (atest +. dtest) in
  Log.info "A_before: %d D_before: %d AD_A_after: %d AD_D_after: %d \
            old: %f new: %f EF: %.3f"
    card_act_before card_dec_before
    card_act_after card_dec_after
    old_rate new_rate (new_rate /. old_rate)

let main () =
  Log.color_on ();
  Log.set_log_level Log.DEBUG;
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n\
              --dataset <string>: path to folder with encoded set molecules\n  \
              --distance <int>: best distance\n  \
              --actives <file>: file with encoded active molecules\n  \
              [--NxCV <int>]: number of folds of cross validation\n  \
              on the training set (default=3)\n  \
              [-np <int>]: number of processors\n  \
              [--dscan <file>]: where to store the scan\n"
        Sys.argv.(0);
      exit 1);
  let ncores = CLI.get_int_def ["-np"] args 1 in
  let dataset_folder = CLI.get_string ["--dataset"] args in
  let actives = CLI.get_string ["--actives"] args in
  let best_d = CLI.get_int_def ["--distance"] args in
  let nfolds = CLI.get_int_def ["--NxCV"] args 3 in
  let dscan_fn = CLI.get_string_def ["--dscan"] args "/dev/null" in
  CLI.finalize();

  let files = Array.to_list (Sys.readdir dataset_folder) in
  List.iter (fun file -> process_file dataset_folder actives file ncores best_d nfolds dscan_fn) files;
  
let () = main ()