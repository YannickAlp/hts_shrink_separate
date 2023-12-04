(* Distance-Based Boolean Applicability Domain (DBBAD)
  reference implementation.

  Copyright (C) 2020, Francois Berenger

  Yamanishi laboratory,
  Department of Bioscience and Bioinformatics,
  Faculty of Computer Science and Systems Engineering,
  Kyushu Institute of Technology,
  680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

  open Printf

  module CLI = Minicli.CLI
  module A = MyArray
  module FpMol = Molenc.FpMol
  module Ht = BatHashtbl
  module L = MyList
  module Log = Dolog.Log
  module Mol = Molecules
  
  let apply_DBBAD ncores best_d chunk active =
    let actives = L.filter FpMol.is_active active in
    let actives_bst = Bstree.(create 10 Two_bands (A.of_list actives)) in
    let maybe_ok_chunk_mols =
      Parany.Parmap.parmap ncores (fun chunk_mol ->
          (Dbad_common.mol_is_inside_global_AD chunk_mol best_d actives_bst,
          chunk_mol)
        ) chunk in
    let ok_chunk_mols =
      L.fold (fun acc (maybe_ok, mol) ->
          if maybe_ok then mol :: acc else acc
        ) [] maybe_ok_chunk_mols in
    let ok_card = L.length ok_chunk_mols in
    let chunk_card = L.length chunk in
    Log.info "passed AD: %d / %d" ok_card chunk_card;
    ok_chunk_mols
  
  let main () =
    Log.color_on ();
    Log.set_log_level Log.DEBUG;
    let argc, args = CLI.init () in
    if argc = 1 then
      (eprintf "usage:\n\
                %s\n\
                --chunk <file>: file with encoded chunk set molecules\n  \
                --active <file>: file with encoded active set molecules\n  \
                [-d;--distance]: best distance
                [-np <int>]: number of processors\n"
        Sys.argv.(0);
      exit 1);
    (*let seed = match CLI.get_int_opt ["--seed"] args with
      | None -> (BatRandom.self_init ();
                BatRandom.int (int_of_float ((2. ** 30.) -. 1.)))
      | Some n -> n in*)
    let ncores = CLI.get_int_def ["-np"] args 1 in
    let best_d = CLI.get_float ["-d";"--distance"] args in
    let chunk_fn = CLI.get_string ["--chunk"] args in
    let active_fn = CLI.get_string ["--active"] args in
    (*let _train_fn = CLI.get_string ["--train"] args in*)
    (*let nfolds = CLI.get_int_def ["--NxCV"] args 3 in*)
    let chunk_dbbad_fn = chunk_fn ^ ".dbbad" in
    (*let dscan_fn = CLI.get_string_def ["--dscan"] args "/dev/null" in*)
    CLI.finalize();
    (* train the DBBAD on the training set *)
     
    let chunk = Molecules.from_file chunk_fn in
    let active = Molecules.from_file active_fn in
    let chunk_DBBAD = apply_DBBAD ncores best_d chunk active in
    Utls.with_out_file chunk_dbbad_fn (fun out ->
        L.iter (FpMol.to_out out) chunk_DBBAD
      );
    Log.info "chunk_DBBAD written to: %s" chunk_dbbad_fn;
    (* report actives proportion before/after DBBAD *)
    let card_act_before, card_dec_before =
      L.filter_counts FpMol.is_active chunk in
    let card_act_after, card_dec_after =
      L.filter_counts FpMol.is_active chunk_DBBAD in
    let old_rate =
      let abefore = float card_act_before in
      let dbefore = float card_dec_before in
      abefore /. (abefore +. dbefore) in
    let new_rate =
      let aafter = float card_act_after in
      let dafter = float card_dec_after in
      aafter /. (aafter +. dafter) in
    Log.info "A_before: %d D_before: %d AD_A_after: %d AD_D_after: %d \
              old: %f new: %f EF: %.3f"
      card_act_before card_dec_before
      card_act_after card_dec_after
      old_rate new_rate (new_rate /. old_rate)
   
   let () = main ()
   