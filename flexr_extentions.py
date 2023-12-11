import os
import pandas as pd
import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GObject
from gi.repository import Gio
from gi.repository import GLib
import coot
import coot_gui_api
import coot_utils
#import coot_gui
from coot_gui import add_simple_action_to_menu, attach_module_menu_button,molecule_chooser_gui_generic,make_store_for_model_molecule_combobox,make_store_for_map_molecule_combobox


#---------------------------------------------------------------------
#  FLEXR interface
#---------------------------------------------------------------------

def add_module_flexr():
    menu = attach_module_menu_button('FLEXR')
    add_simple_action_to_menu(
    menu, "Open","run_flexr_gui",lambda _simple_action, _arg: add_flexr_gui())

def run_flexr(imol,mol,imap,mtz,ringerfile,branching,densitythreshold,singleconf,ringerplots,altlimit):

    from flexr import main
    from src.building import flexr_build
    from src.building import flexr_analysis

    print('Initiating FLEXR...')

    from src.flexrpkg.top_level import args

    ARGS = args()


    if (ringerfile == '/') & (imap is not None):
        try:
            print('Running Ringer...')
            print(mol,mtz)

            # use phenix to remove alt confs
            os.system('phenix.pdbtools %s remove_alt_confs=True' % (mol))
            mol = mol.split('/')[-1][:-4]+"_modified.pdb"
            imolnoconf = coot.read_pdb(mol)
            ARGS.cootmolnum = imolnoconf

            os.system('phenix.maps %s %s' % (mol,mtz))
            map_coef = mol.split()[0][:-4]+"_map_coeffs.mtz"
            print(map_coef)

            os.system('mmtbx.ringer %s %s' % (mol,map_coef))
            ringerfile = mol.split('/')[-1][:-4]+"_ringer.csv"
            ARGS.filename = ringerfile
            altsfile = ringerfile[:-4]+"_"+str(densitythreshold)+"_alts.csv"
        except:
            print('Cannot find Phenix/CCTBX or it failed.')
            print('Done.')
            raise ValueError

    else:
        ARGS.filename = ringerfile
        altsfile = ringerfile[:-4]+"_"+str(densitythreshold)+'_alts.csv'

    # options in the gui
    ARGS.sigmathreshold = densitythreshold
    ARGS.build_limit = altlimit
    ARGS.branching = branching
    ARGS.singleconfs = singleconf
    ARGS.plot = ringerplots

    # these opts so we can call the build function instead of running the script.
    ARGS.pdb = 'None'
    ARGS.build = False

    # run the main flexr script
    main(ARGS)

    # run the building step
    flexrmolnum = flexr_build.building_run(altsfile,ARGS.pdb,ARGS.branching,ARGS.cootmolnum,ARGS.exitcoot)
    #try:
    flexr_analysis.output_summaries(altsfile,imol,flexrmolnum)
    #except:
    #    print('Error producing output summary.')

    print('FLEXR is finished.')
    print('')


def add_flexr_gui():

    def delete_event(*args):
        window.destroy()
        return False

    # get current MODEL for building
    def get_molecule():
        tree_iter = combobox_molecule.get_active_iter()
        imol = -1
        if tree_iter is not None:
            model = combobox_molecule.get_model()
            it = model[tree_iter]
            imol = it[0]
            name = coot.molecule_name(imol)
        return imol,name

    # get current MAP for building
    def get_mtz():
        tree_iter = combobox_map.get_active_iter()
        imap = -1
        if tree_iter is not None:
            mtz = combobox_map.get_model()
            it = mtz[tree_iter]
            imap = it[0]
            name = coot.molecule_name(imap)
        return imap,name

    # some user input, text
    def get_ringer_path():
        t = dist_entry.get_text()
        try:
            ret = str(t)
        except:
            ret = False
        return ret

    # some user input, text
    def get_threshold():
        t = dist_entry2.get_text()
        try:
            ret = float(t)
        except:
            ret = False
        return ret

    # some user input, text
    alt_limit_list = list(range(2,11))
    def get_alt_limit():
        t = alt_limit.get_active()
        try:
            ret = int(t)
            ret = alt_limit_list[ret]
        except:
            ret = False
        return ret

    # Sidechain branching
    # some user input, drop down box
    def get_getbranching():
        at = combobox_coordination.get_active()
        n = str(at)
        return n

    # Single conformer modelling
    # some user input, drop down box
    def get_singleconf():
        at = combobox_coordination2.get_active()
        n = bool(at)
        return n

    # Ringer plots
    # some user input, drop down box
    def get_ringerplots():
        at = combobox_coordination3.get_active()
        n = bool(at)
        return n

    def apply_cb(*args):
        imol,mol = get_molecule()
        imap,mtz = get_mtz()
        ringerpath = get_ringer_path()
        densitythreshold = get_threshold()
        altlimit = get_alt_limit()
        branching = get_getbranching()
        singleconf = get_singleconf()
        ringerplots = get_ringerplots()
        if ringerpath:
            run_flexr(imol,mol,imap,mtz,ringerpath,branching,densitythreshold,singleconf,ringerplots,altlimit)
            #print("now call update_water_results() with imol", imol, "coordination number", n, "imol", imol)
            #update_water_results(imol, n, d)

    # design GUI

    window = Gtk.Window()
    window.set_title("FLEXR: multi-conformer modelling")
    vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    h_sep = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)

    ## header
    hint_text = Gtk.Label()
    hint_text.set_markup("\n For full instructions, go to our <a href=\"https://github.com/TheFischerLab/FLEXR-GUI\" "
                 "title=\"\">GitHub</a> repo")
    vbox.append(hint_text)

    hint_text2 = Gtk.Label()
    hint_text2.set_markup("\n Please cite: <a href=\"https://doi.org/10.1107/S2059798323002498\" "
                "title=\"\">Stachowski and Fischer, Acta Cryst. (2023) D79, 354-357</a> \n")
    vbox.append(hint_text2)

    ## drop down for picking molecule
    hint_text = Gtk.Label(label="\n Molecule:")
    combobox_molecule = Gtk.ComboBox()
    vbox.append(hint_text)
    combobox_mol_items = make_store_for_model_molecule_combobox()
    combobox_molecule.set_model(combobox_mol_items)

    renderer_text = Gtk.CellRendererText()
    if len(combobox_mol_items) > 0:
        combobox_molecule.set_active(0)
    combobox_molecule.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox_molecule.pack_start(renderer_text, True)
    combobox_molecule.add_attribute(renderer_text, "text", 1)
    combobox_molecule.set_margin_start(5)
    combobox_molecule.set_margin_end(5)
    combobox_molecule.set_margin_top(10)
    combobox_molecule.set_margin_bottom(10)
    vbox.append(combobox_molecule)

    # some debugging stuff
    print("debug:: add_flexr_gui(): combobox_mol_items:", combobox_mol_items)
    print("debug:: add_flexr_gui(): combobox_molecule:",  combobox_molecule)

    combobox_molecule.set_active(0)

    ## drop down for picking map
    hint_text = Gtk.Label(label="\n MTZ:")
    combobox_map = Gtk.ComboBox()
    vbox.append(hint_text)
    combobox_mol_items2 = make_store_for_map_molecule_combobox()
    combobox_map.set_model(combobox_mol_items2)

    renderer_text = Gtk.CellRendererText()
    if len(combobox_mol_items2) > 0:
        combobox_map.set_active(0)
    combobox_map.set_entry_text_column(1) # Sets the model column which combo_box
                                      # should use to get strings from to be text_column
    combobox_map.pack_start(renderer_text, True)
    combobox_map.add_attribute(renderer_text, "text", 1)
    combobox_map.set_margin_start(5)
    combobox_map.set_margin_end(5)
    combobox_map.set_margin_top(10)
    combobox_map.set_margin_bottom(10)
    vbox.append(combobox_map)

    # some debugging stuff
    print("debug:: add_flexr_gui(): combobox_mol_items:", combobox_mol_items)
    print("debug:: add_flexr_gui(): combobox_molecule:",  combobox_molecule)

    combobox_map.set_active(0)


    results_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    def clear_previous_results():
        for this_vbox in [results_vbox, metal_results_vbox]:
            child = this_vbox.get_first_child()
            while child is not None:
                next_child = child.get_next_sibling()
                this_vbox.remove(child)
                child = next_child

    window.set_child(vbox)

    # default values for text entry
    hbox_max_dist = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    #dist_label = Gtk.Label(label="FLEXR alt file location: ")
    dist_label = Gtk.Label(label="Ringer CSV (if precomputed): ")
    dist_entry = Gtk.Entry()
    hbox_max_dist.append(dist_label)
    hbox_max_dist.append(dist_entry)
    hbox_max_dist.set_margin_start(6)
    hbox_max_dist.set_margin_end(6)
    hbox_max_dist.set_margin_top(4)
    hbox_max_dist.set_margin_bottom(4)
    vbox.append(hbox_max_dist)
    dist_entry.set_text("/")
    #dist_entry.set_width_chars(5)

    # default values for text entry
    density_threshold = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    density_threshold_label = Gtk.Label(label="Density threshold: ")
    dist_entry2 = Gtk.Entry()
    density_threshold.append(density_threshold_label)
    density_threshold.append(dist_entry2)
    density_threshold.set_margin_start(6)
    density_threshold.set_margin_end(6)
    density_threshold.set_margin_top(4)
    density_threshold.set_margin_bottom(4)
    vbox.append(density_threshold)
    dist_entry2.set_text("0.30")
    #dist_entry.set_width_chars(5)

    # default values for text entry
    #alt_limit = Gtk.ComboBoxText(orientation=Gtk.Orientation.HORIZONTAL)
    alt_limit = Gtk.ComboBoxText.new()
    alt_limit_label = Gtk.Label(label="  Max. number of alts/residue: ")
    #dist_entry3 = Gtk.Entry()
    #alt_limit.append(alt_limit_label)
    #alt_limit.append(dist_entry3)
    for i in alt_limit_list:
        alt_limit.append_text(str(i))
    alt_limit.set_active(1)
    alt_limit.set_margin_start(6)
    alt_limit.set_margin_end(1)
    alt_limit.set_margin_top(4)
    alt_limit.set_margin_bottom(4)
    hbox_chooser = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_number_chooser = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    vbox.append(hbox_chooser)
    hbox_number_chooser.append(alt_limit_label)
    hbox_number_chooser.append(alt_limit)
    vbox.append(hbox_number_chooser)
    #dist_entry3.set_text("3")
    #dist_entry.set_width_chars(5)

    # Create places for drop down option
    # coordination number combobox
    number_text = Gtk.Label(label="  Sidechain branching: ")
    combobox_coordination = Gtk.ComboBoxText.new()
    for i in ['ALL','CA']:
        combobox_coordination.append_text(str(i))

    combobox_coordination.set_active(0)
    combobox_coordination.set_margin_start(6)
    combobox_coordination.set_margin_end(6)
    combobox_coordination.set_margin_top(4)
    combobox_coordination.set_margin_bottom(4)
    hbox_chooser = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_number_chooser = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    vbox.append(hbox_chooser)
    hbox_number_chooser.append(number_text)
    hbox_number_chooser.append(combobox_coordination)
    vbox.append(hbox_number_chooser)

    # Create places for drop down option
    # coordination number combobox
    number_text2 = Gtk.Label(label="  Single conformer modelling (beta): ")
    combobox_coordination2 = Gtk.ComboBoxText.new()
    for i in ['False','True']:
        combobox_coordination2.append_text(str(i))

    combobox_coordination2.set_active(0)
    combobox_coordination2.set_margin_start(6)
    combobox_coordination2.set_margin_end(6)
    combobox_coordination2.set_margin_top(4)
    combobox_coordination2.set_margin_bottom(4)
    hbox_chooser2 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_number_chooser2 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    vbox.append(hbox_chooser2)
    hbox_number_chooser2.append(number_text2)
    hbox_number_chooser2.append(combobox_coordination2)
    vbox.append(hbox_number_chooser2)

    # Create places for drop down option
    # coordination number combobox
    number_text3 = Gtk.Label(label="  Ringer plotting (slow): ")
    combobox_coordination3 = Gtk.ComboBoxText.new()
    for i in ['False','True']:
        combobox_coordination3.append_text(str(i))

    combobox_coordination3.set_active(0)
    combobox_coordination3.set_margin_start(6)
    combobox_coordination3.set_margin_end(6)
    combobox_coordination3.set_margin_top(4)
    combobox_coordination3.set_margin_bottom(4)
    hbox_chooser3 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    hbox_number_chooser3 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
    vbox.append(hbox_chooser3)
    hbox_number_chooser3.append(number_text3)
    hbox_number_chooser3.append(combobox_coordination3)
    vbox.append(hbox_number_chooser3)

    ## Place for 'Output 2'
    #scrolled_win = Gtk.ScrolledWindow()
    #metal_results_scrolled_win = Gtk.ScrolledWindow()
    #metal_results_frame = Gtk.Frame()
    #metal_results_vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    #metal_results_label = Gtk.Label(label="Output 1")
    #metal_results_scrolled_win.set_child(metal_results_frame)
    #metal_results_frame.set_child(metal_results_vbox)
    #vbox.append(metal_results_label)
    #vbox.append(metal_results_scrolled_win)

    ## Place for 'Output 1'
    #water_results_label = Gtk.Label(label="Output 2")
    #scrolled_win.set_child(results_vbox)
    #vbox.append(water_results_label)
    #vbox.append(scrolled_win)
    #vbox.append(h_sep)

    ## Place exit/run buttons
    apply_button = Gtk.Button(label="  Run  ")
    close_button = Gtk.Button(label="  Close  ")
    hbox_buttons = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)

    # create buttons to exit/run
    apply_button.set_margin_start(10)
    apply_button.set_margin_end(10)
    apply_button.set_margin_top(4)
    apply_button.set_margin_bottom(4)
    close_button.set_margin_start(6)
    close_button.set_margin_end(6)
    close_button.set_margin_top(4)
    close_button.set_margin_bottom(4)
    hbox_buttons.set_margin_start(6)
    hbox_buttons.set_margin_end(6)
    hbox_buttons.set_margin_top(4)
    hbox_buttons.set_margin_bottom(4)
    hbox_buttons.append(close_button)
    hbox_buttons.append(apply_button)
    vbox.append(hbox_buttons)
    close_button.connect("clicked", delete_event)
    #dist_entry.connect("activate", entry_activate_event)
    apply_button.connect("clicked", apply_cb)

    window.show()

add_module_flexr()

