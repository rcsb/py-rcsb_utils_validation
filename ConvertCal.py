import xml.etree.ElementTree as ET
from mmcif.api.DataCategory import DataCategory
from mmcif.api.PdbxContainers import DataContainer
from mmcif.io.PdbxReader import PdbxReader
from mmcif.io.PdbxWriter import PdbxWriter
from mmcif.io.IoAdapterCore import IoAdapterCore

import sys
import os

class ConvertXML(object):
    def __init__(self):
        self._tree = None
        self._curContainer = None
        self._pypc = None
        self._model = None
        self._type_map = {}
        self._type_id = 0
        self._resparent = {}
        pass

    def readXml(self, fname):
        self._tree = ET.parse(fname)
        self._root = self._tree.getroot()

    def convertMmCif(self):
        # Get ids
        entry = self._root.find("Entry")
        attr = entry.attrib

        entries = []
        block_id=""
        if "pdbid" in attr:
            entries.append({"id":"PDB", "code":attr["pdbid"]})
            block_id += attr["pdbid"]
        if "emdb_id" in attr:
            entries.append({"id":"EMDB", "code":attr["emdb_id"]})
            if len(block_id) != 0:
                block_id += "_"
            block_id += attr["emdb_id"]
                            
        if self._curContainer is None:
            self._curContainer = DataContainer("valid_" + block_id)

        self._addentries(entries)
        self._processEntry()
        self._processEntities()
        self._processModel()
        self._processEntity()
        self._processInstance()
        self._processEM()
        self._processNMRRestraints()
        self._processDistanceRestraintsAnalysis()
        self._processDihedralRestraintsAnalysis()                        
        self._convertsoftware()

        #self._curContainer.printIt()

    def writeMmCif(self, fpathout):
        myDataList = [self._curContainer]
        with open(fpathout, "w") as ofh:
                pdbxW = PdbxWriter(ofh)
                pdbxW.write(myDataList)

    def _processEntities(self):
        """Handle asym ids, etc"""
        if self._model is None:
            return

        b0 = self._model[0]

        aCat = DataCategory("pdbx_vrpt_entity")
        aCat.appendAttribute("id")
        aCat.appendAttribute("type")
        aCat.appendAttribute("description")        

        c0 = b0.getObj("entity")
        for idx in range(c0.getRowCount()):
            eid = c0.getValue("id", idx)
            etype = c0.getValue("type", idx)
            edesc = c0.getValue("pdbx_description", idx)
            rd = {"id": eid, "type": etype, "description": edesc}
            aCat.append(rd)
        self._curContainer.append(aCat)
            

        aCat = DataCategory("pdbx_vrpt_asym")
        aCat.appendAttribute("label_asym_id")
        aCat.appendAttribute("entity_id")

        c0 = b0.getObj("struct_asym")
        for idx in range(c0.getRowCount()):
            asym = c0.getValue("id", idx)
            entity = c0.getValue("entity_id", idx)
            rd = {"label_asym_id": asym, "entity_id": entity}
            aCat.append(rd)
        self._curContainer.append(aCat)
            
        
    def readModel(self, pathIn):
        io = IoAdapterCore()
        self._model = io.readFile(pathIn, selectList=["entity", "struct_asym"])
        
    def _addentries(self, ents):
        """Adds dictionary of key and id"""
        aCat = DataCategory("pdbx_vrpt_database")
        aCat.appendAttribute("id")
        aCat.appendAttribute("code")
        for e in ents:
            rd = {"id": e["id"], "code": e["code"]}
            aCat.append(rd)
        self._curContainer.append(aCat)

    def _processEntry(self):
        """Process attributes from Entry section"""
        entry = self._root.find("Entry")
        attr = entry.attrib

        # step 0: pdbx_vrpt_
        aCat = DataCategory("pdbx_vrpt_summary")
        atMap = [["entry_id", "pdbid"],
# drop                 ["emdb_id", "emdb_id"],
                 ["PDB_deposition_date", "PDB-deposition-date"],
                 ["PDB_revision_number", "PDB-revision-number"],
                 ["PDB_revision_date", "PDB-revision-date"],

#                 ["PDB_res_high", "PDB-resolution"],
#                 ["PDB_res_low", "PDB-resolution-low"],
#                 ["auth_R_work", "PDB-R"],
#                 ["auth_Rfree", "PDB-Rfree"],

                 ["RNA_suiteness", "RNAsuiteness"],
#                 ["EMDB_resolution", "EMDB-resolution"],
                 ["protein_DNA_RNA_entities", "protein-DNA-RNA-entities"],
                 ["model_CA_only", "CA_ONLY"],
                 ["EMDB_deposition_date", "EMDB-deposition-date"],
                 ["report_creation_date", "XMLcreationDate"],
                 ["attempted_validation_steps", "attemptedValidationSteps"],
                 ["no_ligands_for_mogul", "no-ligands-for-mogul"],
                 ["no_ligands_for_buster_report", "no-ligands-for-buster-report"],
                 ["ligands_for_buster_report", "ligands-for-buster-report"],
                 ["no_percentile_property", "no-percentile-property"],
                 # ["contour_level_primary_map", "contour_level_primary_map"],
                 # ["atom_inclusion_all_atoms", "atom_inclusion_all_atoms"],
                 # ["atom_inclusion_backbone", "atom_inclusion_backbone"],
                 # ["author_provided_fsc_resolution_by_cutoff_0.143", "author_provided_fsc_resolution_by_cutoff_0.143"],
                 # ["author_provided_fsc_resolution_by_cutoff_0.333", "author_provided_fsc_resolution_by_cutoff_0.133"],
                 # ["author_provided_fsc_resolution_by_cutoff_0.5", "author_provided_fsc_resolution_by_cutoff_0.5"],
                 # ["author_provided_fsc_resolution_by_cutoff_halfbit", "author_provided_fsc_resolution_by_cutoff_halfbit"],
                 # ["author_provided_fsc_resolution_by_cutoff_onebit", "author_provided_fsc_resolution_by_cutoff_onebit"],
                 # ["author_provided_fsc_resolution_by_cutoff_threesigma", "author_provided_fsc_resolution_by_cutoff_threesigma"],
                 # ["calculated_fsc_resolution_by_cutoff_0.143", "calculated_fsc_resolution_by_cutoff_0.143"],
                 # ["calculated_fsc_resolution_by_cutoff_0.333", "calculated_fsc_resolution_by_cutoff_0.133"],
                 # ["calculated_fsc_resolution_by_cutoff_0.5", "calculated_fsc_resolution_by_cutoff_0.5"],
                 # ["calculated_fsc_resolution_by_cutoff_halfbit", "calculated_fsc_resolution_by_cutoff_halfbit"],
                 # ["calculated_fsc_resolution_by_cutoff_onebit", "calculated_fsc_resolution_by_cutoff_onebit"],
                 # ["calculated_fsc_resolution_by_cutoff_threesigma", "calculated_fsc_resolution_by_cutoff_threesigma"],
#                 ["nmr_models_consistency_flag", "nmr_models_consistency_flag"], 
#                 ["cyrange_error", "cyrange_error"],
#                 ["cyrange_version", "cyrange_version"],
#                 ["nmrclust_error", "nmrclust_error"],
#                 ["nmrclust_version", "nmrclust_version"],
#                 ["nmrclust_representative_model", "nmrclust_representative_model"],
#                 ["medoid_model", "medoid_model"],
#                 ["nmrclust_number_of_outliers", "nmrclust_number_of_outliers"],
#                 ["nmrclust_number_of_models", "nmrclust_number_of_models"],
#                 ["nmrclust_number_of_clusters", "nmrclust_number_of_clusters"],
#                 ["cyrange_number_of_domains", "cyrange_number_of_domains"],
#                 ["chemical_shift_completeness", "chemical_shift_completeness"],
#                 ["chemical_shift_completeness_full_length", "chemical_shift_completeness_full_length"],
#                 ["panav_version", "panav_version"],
#                 ["rci_version", "rci_version"],
#                 ["shiftchecker_version", "shiftchecker_version"],
                 ]

        rd = {}
        for a in atMap:
            att = a[0]
            if a[1] != "?":
                val = attr.get(a[1], "?")
            else:
                val = atlookup[att]

                
            if att in ["ligands_for_buster_report", "no_ligands_for_mogul", "no_ligands_for_buster_report"]:
                if val == "yes":
                    val = "Y"
                elif val == "no":
                    val = "N"

            if val != "?":
                aCat.appendAttribute(att)
                rd[att] = val

        aCat.append(rd)
        self._curContainer.append(aCat)

        self._processChemicalShiftLists(entry)
        
        
        # Step 1 - create the percentile list - as id will be needed elsewhere
        pbins = attr["percentilebins"]
        # percentilebins="all,2.35,xray" -> xray all refers to only xray
        # percentilebins="all,em" -> all entries in archive, all em
        aCat = DataCategory("pdbx_vrpt_percentile_list")
        aCat.appendAttribute("id")
        aCat.appendAttribute("range")
        aCat.appendAttribute("exp_method")
        sbins = pbins.split(",")
        exp = sbins[-1:][0]
        print("EXP is %s" % exp)
        if exp == "xray":
            self.__absid = "1"
            self.__relid = "2"
            # Absolute vs all x-ray only
            aCat.append([self.__absid, sbins[0], "x-ray"])
            aCat.append([self.__relid, sbins[1], "x-ray"])
            emeth="x-ray"
            
        elif exp == "em":
            self.__absid = "1"
            self.__relid = "2"
            # Absolute all vs em only
            aCat.append([self.__absid, "all", "pdb"])
            aCat.append([self.__relid, "all", "electron microscopy"])
            emeth="electron microscopy"

        elif exp == "nmr":
            self.__absid = "1"
            self.__relid = "2"
            # Absolute all vs NMR onlye
            aCat.append([self.__absid, "all", "pdb"])
            aCat.append([self.__relid, "all", "nmr"])
            emeth="nmr"
        else:
            print("PBINS", pbins, exp)
            sys.exit(1)

        self._curContainer.append(aCat)

        # Create an exptl -- can you tell neutron, EC?
        aCat = DataCategory("pdbx_vrpt_exptl")
        aCat.appendAttribute("ordinal")
        aCat.appendAttribute("method")
        aCat.append(["1", emeth])
        self._curContainer.append(aCat)


        #
        # pdbx_vrpt_percentile_type.id
        # pdbx_vrpt_percentile_type.type
        #
        # Search for type
        
        
        # =====================================
        # pdbx_vrpt_percentile_conditions
        # =====================================
        aCat = DataCategory("pdbx_vrpt_percentile_conditions")
        aCat.appendAttribute("id")
        aCat.appendAttribute("type_id")
        aCat.appendAttribute("rank")
        aCat.appendAttribute("res_high")
        aCat.appendAttribute("res_low")
        aCat.appendAttribute("number_entries_total")
        aCat.appendAttribute("percentile_list_id")


        # pv is a mapping between type/relative id and conditions_id
        self._pvpc = {}
        pmap = [["all_atom_clashscore", self.__absid, "numPDBids-absolute-percentile-clashscore", "?", "?"],
                ["all_atom_clashscore", self.__relid, "numPDBids-relative-percentile-clashscore", "high-resol-relative-percentile-clashscore", "low-resol-relative-percentile-clashscore"],
                ["Ramachandran_outlier_percent", self.__absid, "numPDBids-absolute-percentile-percent-rama-outliers", "?", "?"],
                ["Ramachandran_outlier_percent", self.__relid, "numPDBids-relative-percentile-percent-rama-outliers", "high-resol-relative-percentile-percent-rama-outliers", "low-resol-relative-percentile-percent-rama-outliers"],
                ["rotamer_outliers_percent", self.__absid, "numPDBids-absolute-percentile-percent-rota-outliers", "?", "?"],
                ["rotamer_outliers_percent", self.__relid, "numPDBids-relative-percentile-percent-rota-outliers", "high-resol-relative-percentile-percent-rota-outliers", "low-resol-relative-percentile-percent-rota-outliers"],
                ["R_value_R_free", self.__absid, "numPDBids-absolute-percentile-DCC_Rfree", "?", "?"],
                ["R_value_R_free", self.__relid, "numPDBids-relative-percentile-DCC_Rfree", "high-resol-relative-percentile-DCC_Rfree", "low-resol-relative-percentile-DCC_Rfree"],
                ["RSRZ_outliers_percent", self.__absid, "numPDBids-absolute-percentile-percent-RSRZ-outliers", "?", "?"],
                ["RSRZ_outliers_percent", self.__relid, "numPDBids-relative-percentile-percent-RSRZ-outliers", "high-resol-relative-percentile-percent-RSRZ-outliers", "low-resol-relative-percentile-percent-RSRZ-outliers"],
                ["RNAsuiteness_percent", self.__absid, "numPDBids-absolute-percentile-RNAsuiteness", "?", "?"],
                ["RNAsuiteness_percent", self.__relid, "numPDBids-relative-percentile-RNAsuiteness", "high-resol-relative-percentile-RNAsuiteness", "low-resol-relative-percentile-RNAsuiteness"]
                ]
        cid = 1
        for p in pmap:
            ptype = p[0]
            if ptype not in self._type_map:
                self._type_id += 1
                self._type_map[ptype] = self._type_id
            ptype_id = self._type_map[ptype]
            plist = p[1]
            num = attr.get(p[2], "?")
            if num == "?":
                continue
            if p[3] == "?":
                res_high = "?"
            else:
                res_high = attr.get(p[3], "?")
            if p[4] == "?":
                res_low = "?"
            else:
                res_low = attr.get(p[4], "?")
            aCat.append([str(cid), "?", "?",  res_high, res_low, num, plist])
            self._pvpc[str(ptype_id) + "_" + plist] = str(cid)
            cid += 1

                
        #self._curContainer.append(aCat)
                        


        vmap = [ ["all_atom_clashscore", "absolute-percentile-clashscore", self.__absid],
                 ["all_atom_clashscore", "relative-percentile-clashscore", self.__relid],
                 ["Ramachandran_outlier_percent", "absolute-percentile-percent-rama-outliers", self.__absid],
                 ["Ramachandran_outlier_percent", "relative-percentile-percent-rama-outliers", self.__relid],
                 ["rotamer_outliers_percent", "absolute-percentile-percent-rota-outliers", self.__absid],
                 ["rotamer_outliers_percent", "relative-percentile-percent-rota-outliers", self.__relid],
                 ["R_value_R_free", "absolute-percentile-DCC_Rfree", self.__absid],
                 ["R_value_R_free", "relative-percentile-DCC_Rfree", self.__relid],
                 ["RSRZ_outliers_percent", "absolute-percentile-percent-RSRZ-outliers", self.__absid],
                 ["RSRZ_outliers_percent", "relative-percentile-percent-RSRZ-outliers", self.__relid],                 
                 ["RNAsuiteness_percent", "absolute-percentile-RNAsuiteness", self.__absid],
                 ["RNAsuiteness_percent", "relative-percentile-RNAsuiteness", self.__relid],                 
                 ]

        for v in vmap:
            ptype = v[0]
            if ptype not in self._type_map:
                self._type_id += 1
                self._type_map[ptype] = self._type_id
            ptype_id = self._type_map[ptype]
            
            rank = attr.get(v[1], "?")
            if rank == "?":
                continue
            lid = v[2]
            # Look up parent
            cid = self._pvpc.get(str(ptype_id) + "_" + lid, "?")
            row = int(cid) - 1
            aCat.setValue(ptype_id, "type_id", row)
            aCat.setValue(rank, "rank", row)
            
        self._curContainer.append(aCat)

        # =====================================
        # pdbx_vrpt_summary_density
        # =====================================
        aCat = DataCategory("pdbx_vrpt_summary_diffraction")

        # List of attributes - handle special later
        atMap = [ ["exp_method", "?"],
                  ["ordinal", "?"],
                  ["Babinet_b", "babinet_b"],
                  ["bulk_solvent_b", "bulk_solvent_b"],
                  ["Wilson_B_estimate", "WilsonBestimate"],
                  ["I_over_sigma", "IoverSigma"],
                  ["num_miller_indices", "numMillerIndices"],
                  ["Babinet_k", "babinet_k"],
                  ["bulk_solvent_k", "bulk_solvent_k"],
                  ["Padilla_Yeates_L_mean", "TwinL"],
                  ["Padilla_Yeates_L2_mean", "TwinL2"],
                  ["DCC_R", "DCC_R"],
                  ["DCC_Rfree", "DCC_Rfree"],
                  ["EDS_R", "EDS_R"],
                  ["EDS_res_high", "EDS_resolution"],
                  ["EDS_res_low", "EDS_resolution_low"], # or use eds?
                  
                  ["Wilson_B_aniso", "WilsonBaniso"],
                  ["data_anisotropy", "DataAnisotropy"],
                  ["trans_NCS_details", "TransNCS"],
                  ["B_factor_type", "B_factor_type"],
                  ["acentric_outliers", "acentric_outliers"],
                  ["centric_outliers", "centric_outliers"],
                  ["data_completeness", "DataCompleteness"],
                  ["number_reflns_R_free", "num-free-reflections"],
                  ["percent_free_reflections", "percent-free-reflections"],
                  ["percent_RSRZ_outliers", "percent-RSRZ-outliers"],
                  
                  ["PDB_resolution_high", "PDB-resolution"],
                  ["PDB_resolution_low", "PDB-resolution-low"],
                  ["PDB_R", "PDB-R"],
                  ["PDB_Rfree", "PDB-Rfree"],
                 ]

        atlookup = {"exp_method": emeth, "ordinal" : 1}

        rd = {}
        for a in atMap:
            att = a[0]
            aCat.appendAttribute(att)

            if a[1] != "?":
                val = attr.get(a[1], "?")
                if val == "NotAvailable":
                    val = "?"
            else:
                val = atlookup[att]
            rd[att] = val

            
        aCat.append(rd)

        if emeth not in ["electron microscopy"]:
            self._curContainer.append(aCat)


        ######### summary_nmr
        # =====================================
        # pdbx_vrpt_summary_nmr
        # =====================================
        aCat = DataCategory("pdbx_vrpt_summary_nmr")

        # List of attributes - handle special later
        atMap = [ ["exp_method", "?"],
                  ["ordinal", "?"],
                  ["nmr_models_consistency_flag", "nmr_models_consistency_flag"], 
                  ["nmrclust_representative_model", "nmrclust_representative_model"],
                  ["medoid_model", "medoid_model"],
                  ["nmrclust_number_of_outliers", "nmrclust_number_of_outliers"],
                  ["nmrclust_number_of_models", "nmrclust_number_of_models"],
                  ["nmrclust_number_of_clusters", "nmrclust_number_of_clusters"],
                  ["cyrange_number_of_domains", "cyrange_number_of_domains"],
                  ["chemical_shift_completeness", "chemical_shift_completeness"],
                  ["chemical_shift_completeness_full_length", "chemical_shift_completeness_full_length"],
                 ]
                  
        atlookup = {"exp_method": emeth, "ordinal" : 1}

        rd = {}
        for a in atMap:
            att = a[0]
            aCat.appendAttribute(att)

            if a[1] != "?":
                val = attr.get(a[1], "?")
            else:
                val = atlookup[att]
            rd[att] = val

            
        aCat.append(rd)

        if emeth in ["nmr"]:
            self._curContainer.append(aCat)
        
        # =====================================
        # pdbx_vrpt_summary_em
        # =====================================
        aCat = DataCategory("pdbx_vrpt_summary_em")

        # List of attributes - handle special later
        atMap = [ ["exp_method", "?"],
                  ["ordinal", "?"],
                  ["contour_level_primary_map", "contour_level_primary_map"],
                  ["atom_inclusion_all_atoms", "atom_inclusion_all_atoms"],
                  ["atom_inclusion_backbone", "atom_inclusion_backbone"],
                  ["author_provided_fsc_resolution_by_cutoff_pt_143", "author_provided_fsc_resolution_by_cutoff_0.143"],
                  ["author_provided_fsc_resolution_by_cutoff_pt_333", "author_provided_fsc_resolution_by_cutoff_0.133"],
                  ["author_provided_fsc_resolution_by_cutoff_pt_5", "author_provided_fsc_resolution_by_cutoff_0.5"],
                  ["author_provided_fsc_resolution_by_cutoff_halfbit", "author_provided_fsc_resolution_by_cutoff_halfbit"],
                  ["author_provided_fsc_resolution_by_cutoff_onebit", "author_provided_fsc_resolution_by_cutoff_onebit"],
                  ["author_provided_fsc_resolution_by_cutoff_threesigma", "author_provided_fsc_resolution_by_cutoff_threesigma"],
                  ["calculated_fsc_resolution_by_cutoff_pt_143", "calculated_fsc_resolution_by_cutoff_0.143"],
                  ["calculated_fsc_resolution_by_cutoff_pt_333", "calculated_fsc_resolution_by_cutoff_0.133"],
                  ["calculated_fsc_resolution_by_cutoff_pt_5", "calculated_fsc_resolution_by_cutoff_0.5"],
                  ["calculated_fsc_resolution_by_cutoff_halfbit", "calculated_fsc_resolution_by_cutoff_halfbit"],
                  ["calculated_fsc_resolution_by_cutoff_onebit", "calculated_fsc_resolution_by_cutoff_onebit"],
                  ["calculated_fsc_resolution_by_cutoff_threesigma", "calculated_fsc_resolution_by_cutoff_threesigma"],
                  ["EMDB_resolution", "EMDB-resolution"],

                 ]
                  
        atlookup = {"exp_method": emeth, "ordinal" : 1}

        rd = {}
        for a in atMap:
            att = a[0]
            aCat.appendAttribute(att)

            if a[1] != "?":
                val = attr.get(a[1], "?")
            else:
                val = atlookup[att]
            rd[att] = val

            
        aCat.append(rd)

        if emeth in ["electron microscopy"]:
            self._curContainer.append(aCat)

            # =====================================
        # pdbx_vrpt_summary_geometry
        # =====================================
        aCat = DataCategory("pdbx_vrpt_summary_geometry")

        # List of attributes - handle special later
        atMap = [["ordinal", "?"],
                 ["percent_ramachandran_outliers", "percent-rama-outliers"],
                 ["clashscore", "clashscore"],
                 ["angles_RMSZ", "angles_rmsz"],
                 ["bonds_RMSZ", "bonds_rmsz"],
                 ["num_angles_RMSZ", "num_angles_rmsz"],
                 ["num_bonds_RMSZ", "num_bonds_rmsz"],
                 ["percent_rotamer_outliers", "percent-rota-outliers"]
                  
                 ]
        atlookup = {"ordinal" : 1}

        rd = {}
        for a in atMap:
            att = a[0]
            aCat.appendAttribute(att)

            if a[1] == "?":
                val = atlookup[att]
            else:
                val = attr.get(a[1], '?')
            rd[att] = val


        aCat.append(rd)
        self._curContainer.append(aCat)
        


    def _processModel(self):
        models = self._root.findall("Model")
        if len(models) == 0:
            return

        # ======================================
        # pdbx_vrpt_percentile_entity_view
        # ======================================
        aCat = DataCategory("pdbx_vrpt_model_list")
        aCat.appendAttribute("PDB_model_num")
        aCat.appendAttribute("nmrclust_cluster_id")
        aCat.appendAttribute("nmrclust_representative")
            

        mapping = {"model" : "PDB_model_num",
                   "nmrclust_cluster_id" : "nmrclust_cluster_id",
                   "nmrclust_representative" : "nmrclust_representative"
                   }
        for m in models:
            attr = m.attrib

            rd = {}
            for (atid, val)  in attr.items():
                rd[mapping[atid]] = val
            aCat.append(rd)

        self._curContainer.append(aCat)

            
        
    def _processEntity(self):
        modelledEntityInstances = self._root.findall("ModelledEntityInstance")

        # ======================================
        # pdbx_vrpt_percentile_entity_view
        # ======================================
        aCat = DataCategory("pdbx_vrpt_percentile_entity_view")
        aCat.appendAttribute("ordinal")        
        aCat.appendAttribute("conditions_id")
        aCat.appendAttribute("type_id")
        aCat.appendAttribute("label_asym_id")
        aCat.appendAttribute("PDB_model_num")
        aCat.appendAttribute("entity_id")
        aCat.appendAttribute("auth_asym_id")
        aCat.appendAttribute("rank")
        #               type                    abs/rel id   
        mapping = [["RSRZ_outliers_percent", self.__absid, "absolute_RSRZ_percentile"],
                   ["RSRZ_outliers_percent", self.__relid, "relative_RSRZ_percentile"],
                   ["Ramachandran_outlier_percent", self.__absid, "absolute_rama_percentile"],
                   ["Ramachandran_outlier_percent", self.__relid, "relative_rama_percentile"],
                   ["rotamer_outliers_percent", self.__absid, "absolute_sidechain_percentile"],
                   ["rotamer_outliers_percent", self.__relid, "relative_sidechain_percentile"],
                   ]


        ord = 0
        for mei in modelledEntityInstances:
            attr = mei.attrib
            # print(attr)
            label_asym = attr["said"]
            model = attr["model"]
            entity = attr["ent"]
            auth_asym = attr["chain"]

            for m in mapping:
                ptype = m[0]
                lid = m[1]
                lookup = m[2]
                if ptype not in self._type_map:
                    self._type_id += 1
                    self._type_map[ptype] = self._type_id
                ptype_id = self._type_map[ptype]

                # nonpolymers do not have Ramachandran data
                if lookup not in attr:
                    continue
                
                cid = self._pvpc.get(str(ptype_id) + "_" + lid, "?")
                ord += 1
                rd = { "ordinal":str(ord), "conditions_id" : cid, "type_id": ptype_id, "label_asym_id" : label_asym, "PDB_model_num" : model, "entity_id": entity, "auth_asym_id" : auth_asym}
                rd["rank"] = attr[lookup]
                aCat.append(rd)


        self._curContainer.append(aCat)

        #
        # pdbx_vrpt_percentile_type.id
        # pdbx_vrpt_percentile_type.type
        #
        aCat = DataCategory("pdbx_vrpt_percentile_type")
        aCat.appendAttribute("id")
        aCat.appendAttribute("type")        

        idlist = []
        rev = {}
        for p in self._type_map:
            k = self._type_map[p]
            idlist.append(k)
            rev[k] = p
        idlist.sort()
        for k in idlist:
            rd = {"id": k, "type": rev[k]}
            aCat.append(rd)
            
        self._curContainer.append(aCat)



        # =======================================
        # pdbx_vrpt_summary_entity_geometry
        # =======================================
        aCat = DataCategory("pdbx_vrpt_summary_entity_geometry")
        #               
        mapping = [["ordinal", "?"],
                   ["PDB_model_num", "model"],
                   ["entity_id", "ent"],
                   ["auth_asym_id", "chain"],
                   ["label_asym_id", "said"],                   
                   ["angles_RMSZ", "angles_rmsz"],
                   ["bonds_RMSZ",  "bonds_rmsz"],
                   ["num_angles_RMSZ", "num_angles_rmsz"],
                   ["num_bonds_RMSZ", "num_bonds_rmsz"],
                   ["average_residue_inclusion", "average_residue_inclusion"],
                   ]

        for m in mapping:
            aCat.appendAttribute(m[0])

        ord = 0
        for mei in modelledEntityInstances:
            attr = mei.attrib
            rd = {}
            for m in mapping:
                if m[0] == "ordinal":
                    ord += 1
                    rd[m[0]] = str(ord)
                else:
                    rd[m[0]] = attr.get(m[1], "?")
            aCat.append(rd)

        self._curContainer.append(aCat)
                

    def reslookup(self, model, auth_chain, auth_seq_id, auth_comp_id, auth_ins_code, altcode):
        """Returns lookup string for mapping"""
        if auth_ins_code == " ":
            auth_ins_code = "?"
        return "-".join([model, auth_chain, auth_seq_id, auth_comp_id, auth_ins_code, altcode])
        
    def _processInstance(self):
        """Converts instance level info"""



            

        #mapping from residue to pdbx_vrpt_model_instance.id
        allreskeys = []
        
        # Create a parent category for a residue -- not model specific
        modelledSubGroups = self._root.findall("ModelledSubgroup")

        # ======================================
        # pdbx_vrpt_model_instance
        # ======================================
        aMICat = DataCategory("pdbx_vrpt_model_instance")
        mapping = [["id", "?"],
                   ["PDB_model_num", "model"],
                   ["entity_id", "ent"],
                   ["label_asym_id", "said"],
                   ["label_seq_id", "seq"],
                   ["label_comp_id", "resname"],
                   ["auth_asym_id", "chain"],
                   ["auth_seq_id", "resnum"],
                   ["label_alt_id", "altcode"],
                   ["PDB_ins_code", "icode"]
                   ]

        for m in mapping:
            aMICat.appendAttribute(m[0])

        ord = 0
        for msg in modelledSubGroups:
            attr = msg.attrib
            ord = ord + 1
            rd = {}
            for m in mapping:
                key = m[0]
                lookup = m[1]
                if key == "id":
                    continue
                val = attr[lookup]
                # insertion code
                if val == " ":
                    val = "?"
                rd[key] = val

            altcode = attr["altcode"]
            if altcode == " ":
                altcode = "?"
            reskey = self.reslookup(attr["model"], attr["chain"], attr["resnum"], attr["resname"], attr["icode"], altcode)
            self._resparent[reskey] = str(ord)
            allreskeys.append(str(ord))
            rd["id"] = str(ord)
            aMICat.append(rd)
            
#        self._curContainer.append(aMICat)

        # =====================================
        # pdbx_vrpt_instance_density
        # =====================================
        aCat = DataCategory("pdbx_vrpt_model_instance_density")
        mapping = [["ordinal", "?"],
                   ["instance_id", "?"],
                   ["label_alt_id", "altcode"],
                   ["natoms_eds", "NatomsEDS"],
                   ["RSRCC", "rscc"],
                   ["RSR", "rsr"],                   ["RSRZ", "rsrz"],
                   ["lig_RSRZ_nbr_id", "lig_rsrz_nbr_id"],
                   ["lig_RSR_nbr_mean", "ligRSRnbrMean"],
                   ["lig_RSR_nbr_stdev", "ligRSRnbrStdev"],
                   ["lig_RSR_numnbrs", "ligRSRnumnbrs"],
                   ["lig_RSRZ", "ligRSRZ"],
                   ]

        for m in mapping:
            aCat.appendAttribute(m[0])

        ord = 0
        keep = False

        rsrkeys = {}
        for msg in modelledSubGroups:
            attr = msg.attrib
            # Skip if no rsrz
            if "rsrz" not in attr:
                continue
            ord = ord + 1
            rd = {}
            for m in mapping:
                key = m[0]
                lookup = m[1]
                if lookup == "?":
                    if key == "ordinal":
                        val = str(ord)
                    elif key=="instance_id":
                        altcode = attr["altcode"]
                        if altcode == " ":
                            altcode = "?"
                        reskey = self.reslookup(attr["model"], attr["chain"], attr["resnum"], attr["resname"], attr["icode"], altcode)
                        val = self._resparent[reskey]
                        rsrkeys[val] = str(ord)
                    else:
                        print("UNKNOWN KEY %s" % key)
                        sys.exit(1)
                else:
                    val = attr.get(lookup, "?")
                # empty values
                if val == " ":
                    val = "?"
                rd[key] = val
                keep = True
                
            aCat.append(rd)

        if keep:
            self._trimContainer(aCat)
            self._curContainer.append(aCat)


        # =====================================
        # pdbx_vrpt_instance_geometry
        # =====================================
        aCat = DataCategory("pdbx_vrpt_model_instance_geometry")
        mapping = [["ordinal", "?"],
                   ["instance_id", "?"],
                   ["label_alt_id", "altcode"],
                   ["OWAB", "owab"],
                   ["average_occupancy", "avgoccu"],
                   ["rotamer_class", "rota"],
                   ["phi", "phi"],
                   ["psi", "psi"],
                   ["ramachandran_class", "rama"],
                   ["flippable_sidechain", "flippable-sidechain"],
                   ["RNA_score", "RNAscore"],
                   ["RNA_suite", "RNAsuite"],
                   ["RNA_pucker", "RNApucker"],
#                   ["ligand_num_clashes", "ligand_num_clashes"],
#                   ["ligand_num_symm_clashes", "ligand_num_symm_clashes"],
#                   ["ligand_clashes_outlier", "ligand_clashes_outlier"],
                   ["ligand_chirality_outlier", "ligand_chirality_outlier"],
                   ["cis_peptide", "cis_peptide"],
                   ["cyrange_domain_id", "cyrange_domain_id"],
                   ["validate", "validate"],
                   ["num_H_reduce", "num-H-reduce"],
                   ["mogul_ignore" , "mogul-ignore"],
                   ["mogul_angles_RMSZ", "mogul_angles_rmsz"],
                   ["mogul_bonds_RMSZ", "mogul_bonds_rmsz"],
                   ["mogul_RMSZ_num_angles", "mogul_rmsz_numangles"],
                   ["mogul_RMSZ_num_bonds", "mogul_rmsz_numbonds"],
#                   ["ligand_geometry_outlier", "ligand_geometry_outlier"], 
                   ["ligand_density_outlier", "ligand_density_outlier"],
                   ["residue_inclusion", "residue_inclusion"],
                   ]
        

        for m in mapping:
            aCat.appendAttribute(m[0])

        ord = 0
        keep = False
        geomkeys = {}
        for msg in modelledSubGroups:
            attr = msg.attrib
            # Skip if no occupancy mean
            if "avgoccu" not in attr and "RNAscore" not in attr and "phi" not in attr and "rama" not in attr:
                continue
            ord = ord + 1
            rd = {}
            for m in mapping:
                key = m[0]
                lookup = m[1]
                if lookup == "?":
                    if key == "ordinal":
                        val = str(ord)
                    elif key=="instance_id":
                        altcode = attr["altcode"]
                        if altcode == " ":
                            altcode = "?"
                        reskey = self.reslookup(attr["model"], attr["chain"], attr["resnum"], attr["resname"], attr["icode"], altcode)
                        val = self._resparent[reskey]
                        geomkeys[val] = str(ord)
                    else:
                        print("UNKNOWN KEY %s" % key)
                        sys.exit(1)
                else:
                    val = attr.get(lookup, "?")
                    if lookup in ["ligand_chirality_outlier", "cis_peptide"]:
                        if val == "yes":
                            val = "Y"
                        elif val == "no":
                            val = "N"
#                    if val == "NotAvailable":
#                        val = "?"


                # empty values
                if val == " ":
                    val = "?"
                rd[key] = val
                keep = True

            aCat.append(rd)

        
        if keep:
            # Trim empty attribute columns
            self._trimContainer(aCat)
            self._curContainer.append(aCat)


        def instance_subcat(aCat, mapping, name, splitatoms=False):
            for m in mapping:
                aCat.appendAttribute(m[0])

            ord = 0
            keep = False
            outliervar = {}
            for msg in modelledSubGroups:
                pattr = msg.attrib
                altcode = pattr["altcode"]
                if altcode == " ":
                    altcode = "?"
                reskey = self.reslookup(pattr["model"], pattr["chain"], pattr["resnum"], pattr["resname"], pattr["icode"], altcode)

                aoutlier = msg.findall(name)
                for ao in aoutlier:
                    attr = ao.attrib
                    ord = ord + 1
                    rd = {}
                    for m in mapping:
                        key = m[0]
                        lookup = m[1]
                        if lookup == "?":
                            if key == "ordinal":
                                val = str(ord)
                            elif key=="instance_id":
                                val = self._resparent[reskey]
#                                outliervar[val] = "Y"
                                outliervar[val] = outliervar.get(val, 0) + 1
                            elif key=="label_alt_id":
                                val = pattr["altcode"]
                            elif key=="atom0" and splitatoms:
                                val = attr["atoms"].split(",")[0]
                            elif key=="atom1" and splitatoms:
                                val = attr["atoms"].split(",")[1]
                            elif key=="atom2" and splitatoms:
                                val = attr["atoms"].split(",")[2]
                            elif key=="atom3" and splitatoms:
                                val = attr["atoms"].split(",")[3]
                            elif key=="atom_1" and splitatoms:
                                val = attr["atoms"].split(",")[0]
                            elif key=="atom_2" and splitatoms:
                                val = attr["atoms"].split(",")[1]
                            elif key=="atom_3" and splitatoms:
                                val = attr["atoms"].split(",")[2]
                            elif key=="atom_4" and splitatoms:
                                val = attr["atoms"].split(",")[3]
                            else:
                                print("UNKNOWN KEY %s" % key)
                                sys.exit(1)
                        else:
                            val = attr.get(lookup, "?")
                            # Special case -- change works
                            if lookup == "link":
                                if val == "yes":
                                    val = "Y"
                                elif val == "no":
                                    val = "N"
                                
                        # empty values
                        if val == " ":
                            val = "?"
                        rd[key] = val
                        keep = True

                    aCat.append(rd)

            if keep:
                # aCat.printIt()
                self._curContainer.append(aCat)
                return outliervar
            return {}

        # =====================================
        # pdbx_vrpt_instance_angle_outliers
        # =====================================
        aCat = DataCategory("pdbx_vrpt_instance_intra_angle_outliers")
        mapping = [["ordinal", "?"],
                   ["instance_id", "?"],
                   ["label_alt_id", "?"],
                   ["atom_1", "atom0"],
                   ["atom_2", "atom1"],
                   ["atom_3", "atom2"],
                   ["obs", "obs"],
                   ["mean", "mean"],
                   ["stdev", "stdev"],
                   ["Z", "z"],
                   ["link", "link"],
                   ]

        angleoutliers = instance_subcat(aCat, mapping, "angle-outlier")

        # =====================================
        # pdbx_vrpt_instance_mogul_angle_outliers
        # =====================================
        aCat = DataCategory("pdbx_vrpt_instance_mogul_angle_outliers")
        mapping = [["ordinal", "?"],
                   ["instance_id", "?"],
                   ["label_alt_id", "?"],
#                   ["atoms", "atoms"],
                   ["atom_1", "?"],
                   ["atom_2", "?"],
                   ["atom_3", "?"],                   
                   ["obsval", "obsval"],
                   ["mean", "mean"],
                   ["stdev", "stdev"],
                   ["numobs", "numobs"],                   
                   ["Zscore", "Sscore"],
                   ["mindiff", "mindiff"],
                   ]

        mogangleoutliers = instance_subcat(aCat, mapping, "mog-angle-outlier", splitatoms=True)


        # =====================================
        # pdbx_vrpt_instance_chiral_outliers
        # =====================================
        aCat = DataCategory("pdbx_vrpt_instance_stereo_outliers")
        mapping = [["ordinal", "?"],
                   ["instance_id", "?"],
                   ["label_alt_id", "?"],
                   ["label_atom_id", "atom"],
                   ["problem", "problem"]
                   ]

        chiraloutliers = instance_subcat(aCat, mapping, "chiral-outlier")
            
        # =====================================
        # pdbx_vrpt_instance_bond_outliers
        # =====================================
        aCat = DataCategory("pdbx_vrpt_instance_intra_bond_outliers")
        mapping = [["ordinal", "?"],
                   ["instance_id", "?"],
                   ["label_alt_id", "?"],
                   ["atom_1", "atom0"],
                   ["atom_2", "atom1"],
                   ["obs", "obs"],
                   ["mean", "mean"],
                   ["stdev", "stdev"],
                   ["Z", "z"],
                   ["link", "link"],
                   ]
        

        bondoutliers = instance_subcat(aCat, mapping, "bond-outlier")

        # =====================================
        # pdbx_vrpt_instance_mogul_bond_outliers
        # =====================================
        aCat = DataCategory("pdbx_vrpt_instance_mogul_bond_outliers")
        mapping = [["ordinal", "?"],
                   ["instance_id", "?"],
                   ["label_alt_id", "?"],
                   ["atom_1", "?"],
                   ["atom_2", "?"],                   
                   ["obsval", "obsval"],
                   ["mean", "mean"],
                   ["numobs", "numobs"],
                   ["stdev", "stdev"],
                   ["Zscore", "Zscore"],
                   ["mindiff", "mindiff"],
                   ]

        mogbondoutliers = instance_subcat(aCat, mapping, "mog-bond-outlier", splitatoms=True)

        # =====================================
        # pdbx_vrpt_instance_mogul_torsion_outliers
        # =====================================
        aCat = DataCategory("pdbx_vrpt_instance_mogul_torsion_outliers")
        mapping = [["ordinal", "?"],
                   ["instance_id", "?"],
                   ["label_alt_id", "?"],
                   ["atom_1", "?"],
                   ["atom_2", "?"],                   
                   ["atom_3", "?"],
                   ["atom_4", "?"],                   
                   ["obsval", "obsval"],
                   ["mean", "mean"],
                   ["mindiff", "mindiff"],
                   ["numobs", "numobs"],
                   ["stdev", "stdev"],
                   ["local_density", "local_density"],
                   ]

        mogtorsoutliers = instance_subcat(aCat, mapping, "mog-torsion-outlier",splitatoms=True)


        # =====================================
        # pdbx_vrpt_instance_mogul_ring_outliers
        # =====================================
        aCat = DataCategory("pdbx_vrpt_instance_mogul_ring_outliers")
        mapping = [["ordinal", "?"],
                   ["instance_id", "?"],
                   ["label_alt_id", "?"],
                   ["atoms", "atoms"],
                   ["mean", "mean"],
                   ["mindiff", "mindiff"],
                   ["numobs", "numobs"],
                   ["stdev", "stdev"],
                   ]

        mogringoutliers = instance_subcat(aCat, mapping, "mog-ring-outlier")

        # =====================================
        # pdbx_vrpt_instance_plane_outliers
        # =====================================
        aCat = DataCategory("pdbx_vrpt_instance_intra_plane_outliers")
        mapping = [["ordinal", "?"],
                   ["instance_id", "?"],
                   ["label_alt_id", "?"],
                   ["type", "type"],
                   ["improper", "improper"],
                   ["omega", "omega"],
                   ["plane_rmsd", "planeRMSD"],
                   ]
        

        planeoutliers = instance_subcat(aCat, mapping, "plane-outlier")


        # =====================================
        # pdbx_vrpt_instance_clashes
        # =====================================
        aCat = DataCategory("pdbx_vrpt_instance_clashes")
        mapping = [["ordinal", "?"],
                   ["label_alt_id", "?"],
                   ["instance_id", "?"],
                   ["label_atom_id", "atom"],
                   ["cid", "cid"],
                   ["clashmag", "clashmag"],
                   ["dist", "dist"],
                   ]


        clashoutliers = instance_subcat(aCat, mapping, "clash")

        
        # =====================================
        # pdbx_vrpt_instance_symm_clashes
        # =====================================
        aCat = DataCategory("pdbx_vrpt_instance_symm_clashes")
        mapping = [["ordinal", "?"],
                   ["label_alt_id", "?"],
                   ["instance_id", "?"],
                   ["label_atom_id", "atom"],
                   ["symop", "symop"],
                   ["scid", "scid"],
                   ["clashmag", "clashmag"],
                   ["dist", "dist"],
                   ]

        symclashoutliers = instance_subcat(aCat, mapping, "symm-clash")

        ############## Do we want a lookup index? ######################
        # =====================================
        # pdbx_vrpt_model_instance - add
        # =====================================
#        aCat = DataCategory("pdbx_vrpt_instance_feature")
#        aMICat.appendAttribute("rsr_id")
#        aMICat.appendAttribute("geometry_id")
        aMICat.appendAttribute("count_angle_outliers")
        aMICat.appendAttribute("count_bond_outliers")
        aMICat.appendAttribute("count_clashes")
        aMICat.appendAttribute("count_symm_clashes")
        aMICat.appendAttribute("count_chiral_outliers")
        aMICat.appendAttribute("count_plane_outliers")
        aMICat.appendAttribute("count_mogul_angle_outliers")
        aMICat.appendAttribute("count_mogul_bond_outliers")
        aMICat.appendAttribute("count_mogul_torsion_outliers")
        aMICat.appendAttribute("count_mogul_ring_outliers")        

        for k in allreskeys:
            inst_id = k
            row = int(k) -1
#            aMICat.setValue(rsrkeys.get(k, "?"), "rsr_id", row)
#            aMICat.setValue(geomkeys.get(k, "?"), "geometry_id", row)
            aMICat.setValue(angleoutliers.get(k, "?"), "count_angle_outliers", row)
            aMICat.setValue(bondoutliers.get(k, "?"), "count_bond_outliers", row)
            aMICat.setValue(clashoutliers.get(k, "?"), "count_clashes", row)
            aMICat.setValue(symclashoutliers.get(k, "?"), "count_symm_clashes", row)            
            aMICat.setValue(chiraloutliers.get(k, "?"), "count_chiral_outliers", row)
            aMICat.setValue(planeoutliers.get(k, "?"), "count_plane_outliers", row)
            aMICat.setValue(mogangleoutliers.get(k, "?"), "count_mogul_angle_outliers", row)
            aMICat.setValue(mogbondoutliers.get(k, "?"), "count_mogul_bond_outliers", row)
            aMICat.setValue(mogtorsoutliers.get(k, "?"), "count_mogul_torsion_outliers", row)
            aMICat.setValue(mogringoutliers.get(k, "?"), "count_mogul_ring_outliers", row)


        self._trimContainer(aMICat, nodel=["label_altid", "PDB_ins_code"])

        # Add modelInstance
        self._curContainer.append(aMICat)

    def _convertsoftware(self):
        programs = self._root.find("programs")

        # Create category - we are the only one
        mapping = {"name": "name",
                   "version": "version",
                   "classification": "properties",
                   
                   }

        # details will be keyed based on program name
        keys = ["ordinal", "name", "version", "success_y_or_n", "classification", "details"]
        aCat = DataCategory("pdbx_vrpt_software")
        for k in keys:
            aCat.appendAttribute(k)

        ord = 0
        for program in programs.findall("program"):
            ord = ord+1
            rd = {"ordinal": str(ord),
                  "success_y_or_n" : "?",
                  "details" : "?"}
            
            for att in mapping:
                rd[att] = program.attrib[mapping[att]]
            aCat.append(rd)
            
           #            print(program.attrib)

        # Now add information from Entry
        entry = self._root.find("Entry")
        attr = entry.attrib
        if "cyrange_version" in attr:
            success = "n"
            if attr["cyrange_error"] == "success":
                success = "y"
            ord +=1
            rd = {"ordinal": str(ord),
                  "name" : "cyrange",
                  "version" : attr["cyrange_version"],
                  "classification" : "cyrange",
                  "success_y_or_n" : success
                  }

            aCat.append(rd)

        if "nmrclust_version" in attr:
            success = "n"
            if attr["nmrclust_version"] == "success":
                success = "y"
            ord +=1
            rd = {"ordinal": str(ord),
                  "name" : "nmrclust",
                  "version" : attr["nmrclust_version"],
                  "classification" : "nmrclust",
                  "success_y_or_n" : success
                  }

            aCat.append(rd)


               
        self._curContainer.append(aCat)

    def _processChemicalShiftLists(self, entry):
        css = entry.findall("chemical_shift_list")
        if len(css) == 0:
            return

        # ======================================
        # pdbx_vrpt_chemical_shift_list
        # ======================================
        aCat = DataCategory("pdbx_vrpt_chemical_shift_list")
        aCat.appendAttribute("ordinal")

        alist = ["file_id", "file_name", "block_name", "list_id", "type", "number_of_errors_while_mapping",
                 "number_of_warnings_while_mapping", "number_of_mapped_shifts", "number_of_parsed_shifts",
                 "total_number_of_shifts", "number_of_unparsed_shifts"]

        for a in alist:
            aCat.appendAttribute(a)
            

        ord = 0
        for c in css:
            attr = c.attrib

            ord = ord + 1
            rd = {"ordinal": ord}
            for (atid, val)  in attr.items():
                rd[atid] = val
            aCat.append(rd)

        self._curContainer.append(aCat)


        def create_subtype(aCat, mapping, name):
            aCat.appendAttribute("ordinal")
            for m in mapping:
                aCat.appendAttribute(m[0])


            ord = 0
            keep = False
            for c in css:
                for f in c.findall(name):
                    attr = f.attrib
                    keep = True
                    ord = ord + 1
                    rd = { "ordinal" : str(ord) }
                    for m in mapping:
                        rd[m[0]] = attr[m[1]]
                    aCat.append(rd)

            if keep:
                self._curContainer.append(aCat)
                

        #############################
        # pdbx_vrpt_random_coil_index
        aCat = DataCategory("pdbx_vrpt_random_coil_index")
        mapping =[["auth_asym_id", "chain"],
                  ["rescode", "rescode"],
                  ["auth_seq_id", "resnum"],
                  ["value", "value"]
                  ]
        create_subtype(aCat, mapping, "random_coil_index")

        #############################
        # pdbx_vrpt_unmapped_chemical_shift
        aCat = DataCategory("pdbx_vrpt_unmapped_chemical_shift")
        mapping =[["chain", "chain"],
                  ["rescode", "rescode"],
                  ["resnum", "resnum"],
                  ["atom", "atom"],
                  ["value", "value"],
                  ["error", "error"],
                  ["ambiguity", "ambiguity"],
                  ["disagnostic", "diagnostic"],
                  ]
        create_subtype(aCat, mapping, "unmapped_chemical_shift")

        
        #############################
        # pdbx_vrpt_unparsed_chemical_shift
        aCat = DataCategory("pdbx_vrpt_unparsed_chemical_shift")
        mapping =[["chem_shift_id", "id"],
                  ["chain", "chain"],
                  ["rescode", "rescode"],
                  ["resnum", "resnum"],
                  ["atom", "atom"],
                  ["value", "value"],
                  ["error", "error"],
                  ["ambiguity", "ambiguity"],
                  ["disagnostic", "diagnostic"],
                  ]
        create_subtype(aCat, mapping, "unparsed_chemical_shift")

        #############################
        # pdbx_vrpt_unparsed_chemical_shift
        aCat = DataCategory("pdbx_vrpt_missing_nmrstar_tag")
        mapping =[["nmrstar_tag_description", "nmrstar_tag_description"],
                  ["nmrstar_tag", "nmrstar_tag"]
                  ]
        create_subtype(aCat, mapping, "missing_nmrstar_tag")

        
        #############################
        # pdbx_vrpt_chemical_shift_outlier
        aCat = DataCategory("pdbx_vrpt_chemical_shift_outlier")
        mapping =[["auth_asym_id", "chain"],
                  ["rescode", "rescode"],
                  ["auth_seq_id", "resnum"],
                  ["label_atom_id", "atom"],
                  ["value", "value"],
                  ["zscore", "zscore"],
                  ["prediction", "prediction"]
                  ]
        create_subtype(aCat, mapping, "chemical_shift_outlier")

        #############################
        # pdbx_vrpt_referencing_offset
        aCat = DataCategory("pdbx_vrpt_referencing_offset")
        mapping =[
            ["label_atom_id", "atom"],
            ["uncertainty", "uncertainty"],
            ["precision", "precision"],
            ["value", "value"],
            ["number_of_measurements", "number_of_measurements"]
        ]

        create_subtype(aCat, mapping, "referencing_offset")
        

        #############################
        # pdbx_vrpt_referencing_offset
        aCat = DataCategory("pdbx_vrpt_assign_completeness_well_defined")
        mapping =[
            ["number_of_assigned_shifts", "number_of_assigned_shifts"],
            ["number_of_unassigned_shifts", "number_of_unassigned_shifts"],
            ["number_of_shifts", "number_of_shifts"],
            ["type", "type"],
            ["element", "element"]
        ]

        create_subtype(aCat, mapping, "assignment_completeness_well_defined")
        
        #############################
        # pdbx_vrpt_referencing_offset
        aCat = DataCategory("pdbx_vrpt_assign_completeness_full_length")
        mapping =[
            ["number_of_assigned_shifts", "number_of_assigned_shifts"],
            ["number_of_unassigned_shifts", "number_of_unassigned_shifts"],
            ["number_of_shifts", "number_of_shifts"],
            ["type", "type"],
            ["element", "element"]
        ]

        create_subtype(aCat, mapping, "assignment_completeness_full_length")
        

    def _processEM(self):
        emv = self._root.find("EM_validation")
        if emv is None:
            return

        # EM categories - one or none
        def create_subtype(aCat, mapping, name):
            data =  emv.find(name)
            if data is None:
                return

            aCat.appendAttribute("ordinal")
            for m in mapping:
                aCat.appendAttribute(m[0])


            attr = data.attrib
            ord = 1
            rd = { "ordinal" : str(ord) }
            for m in mapping:
                rd[m[0]] = attr[m[1]]
            aCat.append(rd)

            self._curContainer.append(aCat)

        
            # <xsd:sequence>
            #     <xsd:element ref="rotationally_averaged_power_spectrum" minOccurs="0" maxOccurs="1"/>
            #     <xsd:element ref="volume_estimate" minOccurs="0" maxOccurs="1"/>
            #     <xsd:element ref="atom_inclusion" minOccurs="0" maxOccurs="1"/>
            #     <xsd:element ref="fsc" minOccurs="0" maxOccurs="1"/>
            # </xsd:sequence>

        #############################
        # pdbx_vrpt_em_details
        aCat = DataCategory("pdbx_vrpt_em_details")
        mapping =[
            ["recommended_contour_level", "value"]
        ]

        create_subtype(aCat, mapping, "RecommendedContourLevel")
        

        # Handle curves.
        #_pdbx_vrpt_em_graph_map_value_distribution.graph_id   map_value_distribution 
        #_pdbx_vrpt_em_graph_rotationally_averaged_power_spectrum.graph_id   rotationally_averaged_power_spectrum 
        # _pdbx_vrpt_em_graph_volume_estimate.graph_id   volume_estimate 
        #loop_
        #_pdbx_vrpt_em_graph_atom_inclusion.graph_id 
        #_pdbx_vrpt_em_graph_atom_inclusion.type 
        #atom_inclusion_all_atoms all_atoms 
        #atom_inclusion_backbone  backbone

        # Parse data, prepare structures

        #         loop_
        # _pdbx_vrpt_em_2d_graph_info.graph_data_id 
        # _pdbx_vrpt_em_2d_graph_info.graph_id 
        # _pdbx_vrpt_em_2d_graph_info.title 
        # _pdbx_vrpt_em_2d_graph_info.x_axis_title 
        # _pdbx_vrpt_em_2d_graph_info.x_axis_scale 
        # _pdbx_vrpt_em_2d_graph_info.x_axis_units 
        # _pdbx_vrpt_em_2d_graph_info.y_axis_title 
        # _pdbx_vrpt_em_2d_graph_info.y_axis_scale 
        # _pdbx_vrpt_em_2d_graph_info.y_axis_units 
        twoDGInfo = DataCategory("pdbx_vrpt_em_2d_graph_info")
        twoDGInfo.appendAttribute("graph_data_id")
        twoDGInfo.appendAttribute("graph_id")
        twoDGInfo.appendAttribute("title")
        twoDGInfo.appendAttribute("x_axis_title")
        twoDGInfo.appendAttribute("x_axis_scale")
        twoDGInfo.appendAttribute("x_axis_units")                
        twoDGInfo.appendAttribute("y_axis_title")
        twoDGInfo.appendAttribute("y_axis_scale")
        twoDGInfo.appendAttribute("y_axis_units")                

        # 
        #loop_
        #_pdbx_vrpt_em_2d_graph_data.ordinal 
        #_pdbx_vrpt_em_2d_graph_data.graph_data_id 
        #_pdbx_vrpt_em_2d_graph_data.x_value 
        #_pdbx_vrpt_em_2d_graph_data.y_value 
        # 1    d_mvd   -0.3271170855 1.4771212547  

        self.twoDGord = 0

        twoDG = DataCategory("pdbx_vrpt_em_2d_graph_data")
        twoDG.appendAttribute("ordinal")
        twoDG.appendAttribute("graph_data_id")
        twoDG.appendAttribute("x_value")
        twoDG.appendAttribute("y_value")


        def getMapInfo(graph_data_id, graph_id, attr):
            rd = {"graph_data_id" : graph_data_id, "graph_id" : graph_id}
            mapping = {
                "Title" : "title",
                "xTitle" : "x_axis_title",
                "xScale" : "x_axis_scale",
                "xUnit" : "x_axis_units",                                 
                "yTitle" : "y_axis_title",
                "yScale" : "y_axis_scale",
                "yUnit" : "y_axis_units"
                }
            for (k, v) in mapping.items():
                if k in attr:
                    rd[v] = attr[k]
            print(rd)
            return rd

        def addGraphInfo(mapinfo):
            twoDGInfo.append(mapinfo)

            
        def addGraph(graph_data_id, coords):
            for c in coords:
                attr = c.attrib
                self.twoDGord += 1
                rd = {"ordinal": str(self.twoDGord),
                      "graph_data_id": graph_data_id,
                      "x_value": attr["x"],
                      "y_value": attr["y"]
                      }
                twoDG.append(rd)


        earlyCat = []

        def doGraph(xmlkey, short):
            mv = emv.findall(xmlkey)

            if mv is  None:
                return


            gid = xmlkey

            aCat = DataCategory("pdbx_vrpt_em_graph_" + xmlkey)
            aCat.appendAttribute("graph_id")
            aCat.append([gid])

            ord = 0
            for m in mv:
                if len(mv) > 1:
                    ord += 1
                    gdid = short + "_" + str(ord)
                else:
                    gdid = short
                
                #
                # Get map info
                attr = m.attrib
                mapinfo = getMapInfo(gdid, gid, attr)
                addGraphInfo(mapinfo)

                coords = m.findall("coordinate")
                addGraph(gdid, coords)

            # REMOVED AS SINGLE INSTANCE
            #earlyCat.append(aCat)

        def doGraph_ai(xmlkey, short):
            mv = emv.find(xmlkey)

            if mv is  None:
                return


            gid = xmlkey

            aCat = DataCategory("pdbx_vrpt_em_graph_" + xmlkey)
            aCat.appendAttribute("graph_id")
            aCat.appendAttribute("type")            

            for m in mv:
                print(m.tag)
                tag = m.tag
                lookup = { "all_atoms" : "aa", "backbone" : "bb"}
                gid = xmlkey + "_" + tag
                gdid = short + "_" + lookup[tag]
                rd = { "graph_id" : gid, "type" : tag}
                aCat.append(rd)
                #
                # Get map info
                attr = m.attrib
                mapinfo = getMapInfo(gdid, gid, attr)
                addGraphInfo(mapinfo)

                coords = m.findall("coordinate")
                addGraph(gdid, coords)

            earlyCat.append(aCat)
            

        def doGraph_fsc(xmlkey, short):
            mv = emv.find(xmlkey)

            if mv is  None:
                return

            # Intersection
            intersections = mv.find("resolution_intersections")
            if intersections is not None:
                aCat = DataCategory("pdbx_vrpt_em_resolution_intersections")
                aCat.appendAttribute("ordinal")
                base = ["ordinal", "resolution_units", "spatial_frequency_units"]
                map = [ "correlation", "resolution", "spatial_frequency", "curve", "type"]
                for m in base + map:
                    aCat.appendAttribute(m)

                resunit = intersections.attrib["resolution_unit"]
                spatial_freq = intersections.attrib["spatial_frequency_unit"]

                ord = 0
                for inter in intersections.findall("intersection"):
                    ord = ord + 1
                    rd = {"ordinal" : str(ord), "resolution_units" : resunit, "spatial_frequency_units" : spatial_freq}
                    for m in map:
                        rd[m] = inter.attrib.get(m, "?")
                    print(rd)
                    aCat.append(rd)
                    
                earlyCat.append(aCat)
                
            # fsc curve
            aCat = DataCategory("pdbx_vrpt_em_graph_fsc_curve")
            aCat.appendAttribute("graph_id")
            aCat.appendAttribute("type")
            aCat.appendAttribute("curve_name")            
            curves = mv.find("fsc_curves") # Must be one
            fsc_curves = curves.findall("fsc_curve")

            ord = 0
            for f in fsc_curves:
                attr = f.attrib
                rd = {"graph_id" : attr["curve_name"],
                      "type" : attr["type"],
                      "curve_name" : attr["curve_name"]
                      }
                aCat.append(rd)

                ord += 1
                
                # Get map info
                gid = attr["curve_name"]
                gdid= "fsc_" + str(ord)

                mapinfo = getMapInfo(gdid, gid, attr)
                addGraphInfo(mapinfo)

                coords = f.findall("coordinate")
                addGraph(gdid, coords)
                
            earlyCat.append(aCat)
                

            # fsc_indicator_curves>

            aCat = DataCategory("pdbx_vrpt_em_graph_fsc_indicator_curve")
            aCat.appendAttribute("graph_id")
            aCat.appendAttribute("type")
            aCat.appendAttribute("curve_name")
            aCat.appendAttribute("data_curve_name")
            
            curves = mv.find("fsc_indicator_curves") # Must be one
            fsc_curves = curves.findall("fsc_indicator_curve")

            ord = 0
            for f in fsc_curves:
                attr = f.attrib
                rd = {"graph_id" : attr["curve_name"],
                      "type" : attr["type"],
                      "curve_name" : attr["curve_name"],
                      "data_curve_name" : attr["data_curve"]                      
                      }
                aCat.append(rd)

                ord += 1
                
                # Get map info
                gid = attr["curve_name"]
                gdid= "fsc_i_" + str(ord)

                mapinfo = getMapInfo(gdid, gid, attr)
                addGraphInfo(mapinfo)

                coords = f.findall("coordinate")
                addGraph(gdid, coords)

                
            earlyCat.append(aCat)
                


            return 


            mv = mv.getchildren()

            print("CHILDREN", mv)
            for m in mv:
                tag = m.tag
                lookup = { "all_atoms" : "aa", "backbone" : "bb"}
                gid = xmlkey + "_" + tag
                gdid = short + "_" + lookup[tag]
                rd = { "graph_id" : gid, "type" : tag}
                aCat.append(rd)
                #
                # Get map info
                attr = m.attrib
                mapinfo = getMapInfo(gdid, gid, attr)
                addGraphInfo(mapinfo)

                coords = m.findall("coordinate")
                addGraph(gdid, coords)

            earlyCat.append(aCat)
            

        doGraph("map_value_distribution", "d_dvd")
        doGraph("rotationally_averaged_power_spectrum", "d_raps")
        doGraph("volume_estimate", "d_ve")
        doGraph_ai("atom_inclusion", "d_ai")
        doGraph_fsc("fsc", "fsc")
        

        for e in earlyCat:
            self._curContainer.append(e)
        self._curContainer.append(twoDGInfo)
        self._curContainer.append(twoDG)            

    def _processNMRRestraints(self):
        nra = self._root.find("NMR_restraints_analysis")
        if nra is None:
            return

        crr = nra.find("conformationally_restricting_restraints")
        if crr:
            
            aCat = DataCategory("pdbx_vrpt_restraint_summary")
            aCat.appendAttribute("ordinal")            
            aCat.appendAttribute("description")
            aCat.appendAttribute("value")

            ord = 0
            for c in crr.findall("restraint_summary"):
                ord += 1
                rd = {"ordinal": str(ord), "description": c.get("description"),
                      "value": c.get("value")}
                aCat.append(rd)

            self._curContainer.append(aCat)


        crr = nra.find("residual_distance_violations")
        if crr:
            
            aCat = DataCategory("pdbx_vrpt_residual_distance_violations")
            keys = ["ordinal", "max_violation", "bins", "violations_per_model"]
            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("residual_distance_violation"):
                rd = {}
                for k in keys:
                    if k == "ordinal":
                        ord += 1
                        rd[k] = str(ord)
                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)

        crr = nra.find("residual_angle_violations")
        if crr:
            
            aCat = DataCategory("pdbx_vrpt_residual_angle_violations")
            keys = ["ordinal", "max_violation", "bins", "violations_per_model"]
            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("residual_angle_violation"):
                rd = {}
                for k in keys:
                    if k == "ordinal":
                        ord += 1
                        rd[k] = str(ord)
                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)


    def _processDistanceRestraintsAnalysis(self):
        nra = self._root.find("distance_restraints_analysis")
        if nra is None:
            return

        crr = nra.find("distance_violations_summary")
        if crr:
            
            aCat = DataCategory("pdbx_vrpt_distance_violation_summary")
            keys = ["ordinal", "restraint_type", "restraint_sub_type", "consistently_violated_count", "consistently_violated_percent_total", "consistently_violated_percent_type", "restraints_count",
                    "violated_count", "percent_total", "violated_percent_total", "violated_percent_type"]
            
            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("distance_violation_summary"):
                rd = {}
                for k in keys:
                    if k == "ordinal":
                        ord += 1
                        rd[k] = str(ord)
                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)

        crr = nra.find("distance_violations_in_models")
        if crr:
            aCat = DataCategory("pdbx_vrpt_distance_violation_model_summary")
            keys = ["ordinal", "PDB_model_num", "max_violation", "mean_violation", "median_violation", "standard_deviation"]
            
            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("distance_violations_in_model"):
                rd = {}
                for k in keys:
                    if k == "ordinal":
                        ord += 1
                        rd[k] = str(ord)
                    elif k=="PDB_model_num":
                        rd[k] = c.get("model")
                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)


            # dist_rest_types
            aCat = DataCategory("pdbx_vrpt_distance_violation_model_restraints")
            keys = ["ordinal", "PDB_model_num", "dist_rest_type", "violations_count"]
            
            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("distance_violations_in_model"):
                model = c.get("model")
                print("MODEL %s" % model)
                for m in c.findall("dist_rest_types"):
                    rd = {}
                    for k in keys:
                        if k=="PDB_model_num":
                            rd[k] = model
                        elif k=="ordinal":
                            ord += 1
                            rd[k] = str(ord)
                        else:
                            rd[k] = m.get(k)
                    aCat.append(rd)

            self._curContainer.append(aCat)
            
        crr = nra.find("distance_violations_in_ensemble")
        if crr:
            aCat = DataCategory("pdbx_vrpt_distance_violations_ensemble_summary")
            keys = ["id", "fraction_of_ensemble_count", "fraction_of_ensemble_percent"]
            
            for k in keys:
                aCat.appendAttribute(k)

            ordinal = 0
            for c in crr.findall("distance_violation_in_ensemble"):
                ordinal = ordinal + 1
                rd = {}
                for k in keys:
                    if k == "id":
                        rd[k] = ordinal
                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)


            # dist_rest_types
            aCat = DataCategory("pdbx_vrpt_distance_violation_ensemble")
            keys = ["ordinal", "ensemble_distance_count", "dist_rest_type", "violations_count"]
            
            for k in keys:
                aCat.appendAttribute(k)

            ordinal = 0
            ord = 0
            for c in crr.findall("distance_violation_in_ensemble"):
                ordinal += 1
                for m in c.findall("dist_rest_types"):
                    rd = {}
                    for k in keys:
                        if k=="ensemble_distance_count":
                            rd[k] = ordinal
                        elif k=="ordinal":
                            ord += 1
                            rd[k] = str(ord)
                        else:
                            rd[k] = m.get(k)
                    aCat.append(rd)

            self._curContainer.append(aCat)
            
        crr = nra.find("most_violated_distance_restraints")
        if crr:
            aCat = DataCategory("pdbx_vrpt_most_violated_distance_restraints")
            keys = ["ordinal", "altcode_1", "altcode_2", "atom_1", "atom_2", "chain_1", "chain_2", "ent_1", "ent_2", "icode_1", "icode_2", "mean_distance_violation",  "median_violation", "resname_1", "resname_2", "resnum_1",  "resnum_2", "rest_id",  "rlist_id", "said_1", "said_2", "seq_1", "seq_2", "standard_deviation", "violated_models"]            
            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("most_violated_distance_restraint"):
                rd = {}
                for k in keys:
                    if k == "ordinal":
                        ord += 1
                        rd[k] = str(ord)
                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)

        crr = nra.find("violated_distance_restraints")
        if crr:
            aCat = DataCategory("pdbx_vrpt_violated_distance_restraints")
            keys = ["ordinal",  "atom_1", "atom_2", "rest_id",  "rlist_id", "violation", "rlist_id", "instance_id_1", "instance_id_2"]            
            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("violated_distance_restraint"):
                rd = {}
                for k in keys:
                    if k == "ordinal":
                        ord += 1
                        rd[k] = str(ord)
                    elif k=="PDB_model_num":
                        rd[k] = c.get("model")
                    elif k == "instance_id_1":
                        altcode = c.get("altcode_1")
                        if altcode == " ":
                            altcode = "?"
                        if altcode == ".":
                            # because NMR is different
                            altcode = "?"

                        reskey = self.reslookup(c.get("model"), c.get("chain_1"), c.get("resnum_1"), c.get("resname_1"), c.get("icode_1"), altcode)
                        rd[k] = self._resparent[reskey]
                    elif k == "instance_id_2":
                        altcode = c.get("altcode_2")
                        if altcode == " ":
                            altcode = "?"
                        if altcode == ".":
                            # because NMR is different
                            altcode = "?"

                        reskey = self.reslookup(c.get("model"), c.get("chain_2"), c.get("resnum_2"), c.get("resname_2"), c.get("icode_2"), altcode)
                        rd[k] = self._resparent[reskey]
                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)
            
    def _processDihedralRestraintsAnalysis(self):
        nra = self._root.find("dihedralangle_restraints_analysis")
        if nra is None:
            return

        crr = nra.find("dihedralangle_violations_summary")
        if crr:
            
            aCat = DataCategory("pdbx_vrpt_dihedralangle_violations_summary")
            keys = ["ordinal", "consistently_violated_count", "consistently_violated_percent_total", "consistently_violated_percent_type", "percent_total", "restraint_type", "restraints_count", "violated_count", "violated_percent_total", "violated_percent_type"]

            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("dihedralangle_violation_summary"):
                rd = {}
                for k in keys:
                    if k == "ordinal":
                        ord += 1
                        rd[k] = str(ord)
                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)
            
        crr = nra.find("dihedralangle_violations_in_models")
        if crr:
            aCat = DataCategory("pdbx_vrpt_dihedralangle_violation_model_summary")
            keys = ["ordinal", "PDB_model_num", "max_violation", "mean_violation", "median_violation", "standard_deviation"]
            
            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("dihedralangle_violations_in_model"):
                rd = {}
                for k in keys:
                    if k=="ordinal":
                        ord += 1
                        rd[k] = str(ord)
                    elif k=="PDB_model_num":
                        rd[k] = c.get("model")
                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)


            # dist_rest_types
            aCat = DataCategory("pdbx_vrpt_dihedralangle_violation_model")
            keys = ["ordinal", "PDB_model_num", "ang_rest_type", "violations_count"]
            
            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("dihedralangle_violations_in_model"):
                model = c.get("model")
                for m in c.findall("ang_rest_types"):
                    rd = {}
                    for k in keys:
                        if k in ["model" "PDB_model_num"]:
                            rd[k] = model
                        elif k=="ordinal":
                            ord += 1
                            rd[k] = str(ord)
                        else:
                            rd[k] = m.get(k)
                    aCat.append(rd)

            self._curContainer.append(aCat)
            
        crr = nra.find("dihedralangle_violations_in_ensemble")
        if crr:
            aCat = DataCategory("pdbx_vrpt_dihedralangle_violation_ensemble_summary")
            keys = ["id", "fraction_of_ensemble_count", "fraction_of_ensemble_percent"]
            for k in keys:
                aCat.appendAttribute(k)

            ordinal = 0
            for c in crr.findall("dihedralangle_violation_in_ensemble"):
                ordinal += 1
                rd = {}
                for k in keys:
                    if k == "id":
                        rd[k] = ordinal
                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)


            # dist_rest_types
            aCat = DataCategory("pdbx_vrpt_dihedralangle_ensemble_violation")
            keys = ["ordinal", "ensemble_dihedral_count", "ang_rest_type", "violations_count"]
            
            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            ordinal = 0
            for c in crr.findall("dihedralangle_violation_in_ensemble"):
                ordinal += 1
                for m in c.findall("ang_rest_types"):
                    rd = {}
                    for k in keys:
                        if k=="ensemble_dihedral_count":
                            rd[k] = ordinal
                        elif k=="ordinal":
                            ord += 1
                            rd[k] = ord
                        else:
                            rd[k] = m.get(k)
                    aCat.append(rd)

            self._curContainer.append(aCat)

        crr = nra.find("most_violated_dihedralangle_restraints")
        if crr:
            aCat = DataCategory("pdbx_vrpt_most_violated_dihedral_restraints")
            keys = ["altcode_1", "chain_1", "resnum_1", "resname_1", "icode_1", "ent_1", "said_1", "seq_1", "atom_1", "chain_2", "altcode_2",
                    "resnum_2", "resname_2", "seq_2", "said_2", "ent_2", "icode_2", "atom_2", "MeanAngleViolation", "altcode_3", "chain_3", "resnum_3",
                    "resname_3", "icode_3", "ent_3", "said_3", "seq_3", "atom_3", "altcode_4", "chain_4", "resnum_4", "resname_4", "icode_4",
                    "ent_4", "said_4", "seq_4", "atom_4", "violated_models"
                    ]

            for k in keys:
                aCat.appendAttribute(k)

            for c in crr.findall("most_violated_dihedralangle_restraint"):
                rd = {}
                for k in keys:
                    rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)
            
        crr = nra.find("violated_dihedralangle_restraints")
        if crr:
            aCat = DataCategory("pdbx_vrpt_violated_dihedral_restraints")
            keys = ["ordinal", "rest_id", "rlist_id", "instance_id_1", "instance_id_2", "instance_id_3", "instance_id_4",
                    "atom_1",
                    "atom_2",
                    "atom_3",
                    "atom_4",
                    "violation"
                    ]

            for k in keys:
                aCat.appendAttribute(k)

            ord = 0
            for c in crr.findall("violated_dihedralangle_restraint"):
                rd = {}
                for k in keys:
                    if k == "ordinal":
                        ord += 1
                        rd[k] = str(ord)
                    elif k=="PDB_model_num":
                        rd[k] = c.get("model")
                    elif k == "instance_id_1":
                        altcode = c.get("altcode_1")
                        if altcode == " ":
                            altcode = "?"
                        if altcode == ".":
                            # because NMR is different
                            altcode = "?"

                        reskey = self.reslookup(c.get("model"), c.get("chain_1"), c.get("resnum_1"), c.get("resname_1"), c.get("icode_1"), altcode)
                        rd[k] = self._resparent[reskey]
                    elif k == "instance_id_2":
                        altcode = c.get("altcode_2")
                        if altcode == " ":
                            altcode = "?"
                        if altcode == ".":
                            # because NMR is different
                            altcode = "?"
                        reskey = self.reslookup(c.get("model"), c.get("chain_2"), c.get("resnum_2"), c.get("resname_2"), c.get("icode_2"), altcode)
                        rd[k] = self._resparent[reskey]

                    elif k == "instance_id_3":
                        altcode = c.get("altcode_3")
                        if altcode == " ":
                            altcode = "?"
                        if altcode == ".":
                            # because NMR is different
                            altcode = "?"

                        reskey = self.reslookup(c.get("model"), c.get("chain_3"), c.get("resnum_3"), c.get("resname_3"), c.get("icode_3"), altcode)
                        rd[k] = self._resparent[reskey]
                    elif k == "instance_id_4":
                        altcode = c.get("altcode_4")
                        if altcode == " ":
                            altcode = "?"
                        if altcode == ".":
                            # because NMR is different
                            altcode = "?"
                        reskey = self.reslookup(c.get("model"), c.get("chain_4"), c.get("resnum_4"), c.get("resname_4"), c.get("icode_4"), altcode)
                        rd[k] = self._resparent[reskey]

                    else:
                        rd[k] = c.get(k)
                aCat.append(rd)

            self._curContainer.append(aCat)

            
    def _trimContainer(self, aCat, nodel=None):
        """Remove columns that have not data"""
        if nodel is None:
            nodel = []
        aset = set(aCat.getAttributeList())

        keep = set()
        for row in range(aCat.getRowCount()):
            rd = aCat.getRowAttributeDict(row)
            for (key, item) in rd.items():
                if item != "?":
                    keep.add(key)
        delist = list(aset-keep)
        print("Deleting %s" % delist)
        for d in delist:
            if d not in nodel:
                aCat.removeAttribute(d)

def main():

    #  Generic, trna, ligand, em, nmr
    reps = ["5b1l", "1ivs", "1o08", "5a32", "6um9", "7jia", "D_1292113885", "2juw", "EMD-21926"]

    for r in reps:
        fnamein = "%s_validation.xml" % r
        modelin = "%s.cif" % r        

        cx = ConvertXML()
        cx.readXml(fnamein)
        if os.path.exists(modelin):
            cx.readModel(modelin)
        cx.convertMmCif()
        cx.writeMmCif("%s_val.cif" % r)
if __name__ == "__main__":
    main()
