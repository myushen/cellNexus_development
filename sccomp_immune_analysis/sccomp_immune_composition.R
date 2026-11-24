
job::job({
  
  library(tidytargets)
  library(targets)
  library(crew)
  library(crew.cluster)
  library(glue)
  
  #-----------------------#
  # FUNCTIONS
  #-----------------------#
  
  result_directory <- "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets"
  dir.create(result_directory, showWarnings = FALSE, recursive = TRUE)
  
  check_rclone_installation <- function() {
    rclone_path <- path.expand("/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone")
    if (!file.exists(rclone_path)) {
      stop("rclone was not found in the expected location '/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone'.")
    }
  }
  
  rclone_copy <- function(from, to) {
    check_rclone_installation()
    system(glue("/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone copy {from} {to}"))
  }
  
  rclone_read_csv <- function(remote, dest = tempdir(), ...) {
    rclone_copy(remote, dest)
    read_csv(file.path(dest, basename(remote)), ...)
  }
  
  edit_covariates = function(tbl){
    
    tissue_grouped = list(
      
      # Respiratory System
      "respiratory system" = c(
        "lung", "lung parenchyma", "alveolus of lung",  "bronchus",
        "respiratory airway", "pleura", "pleural effusion", "middle lobe of right lung",
        "upper lobe of left lung", "lower lobe of left lung", "upper lobe of right lung",
        "lower lobe of right lung", "lingula of left lung", "right lung", "left lung"
      ),
      
      "trachea" = c( "epithelium of trachea", "trachea"),
      
      # Cardiovascular System
      "cardiovascular system" = c(
        "heart", "heart left ventricle", "heart right ventricle", "cardiac ventricle",
        "cardiac atrium", "right cardiac atrium", "left cardiac atrium", "apex of heart",
        "aorta", "coronary artery", 
        "venous blood", "anterior wall of left ventricle", "myocardium", "interventricular septum", "ventricular tissue", "basal zone of heart"
      ),
      
      "vasculature" = c("kidney blood vessel", "artery", "vein", "vasculature", "mesenteric artery"),
      
      # Umbilical Cord Blood
      "umbilical cord blood" = "umbilical cord blood",
      
      # Oesophagus
      "oesophagus" = c(
        "esophagus", "lower esophagus", "esophagus muscularis mucosa",
        "submucosal esophageal gland",
        
        # Epithelium
        "epithelium of esophagus"
      ),
      
      # Stomach
      "stomach" = c(
        "stomach", "body of stomach", "cardia of stomach"
      ),
      
      # Small Intestine
      "small intestine" = c(
        "small intestine", "duodenum", "jejunum", "ileum",
        
        # Epithelium
        "epithelium of small intestine", "jejunal epithelium", "ileal epithelium",
        "submucosa of ileum", "lamina propria of small intestine"
      ),
      
      # Large Intestine
      "large intestine" = c(
        "large intestine", "colon", "left colon", "right colon",
        "sigmoid colon", "descending colon", "transverse colon",
        "ascending colon", "hepatic flexure of colon", "caecum",
        "rectum", "appendix", "vermiform appendix",
        
        # epithelium
        "colonic epithelium", "submucosa of ascending colon", "lamina propria of large intestine",
        "mucosa of colon", "lamina propria of mucosa of colon", "caecum epithelium"
      ),
      
      # Digestive System (General)
      "digestive system (general)" = c(
        "intestine", "hindgut", "lamina propria", "mucosa"
      ),
      
      # Nasal, Oral, and Pharyngeal Regions
      "nasal, oral, and pharyngeal regions" = c(
        "nasal cavity", "nasopharynx", "oral mucosa", "tongue", "anterior part of tongue",
        "posterior part of tongue", "gingiva", "nose", "saliva"
      ),
      
      # Cerebral Lobes and Cortical Areas
      "cerebral lobes and cortical areas" = c(
        "frontal lobe", "left frontal lobe", "right frontal lobe", "primary motor cortex",
        "dorsolateral prefrontal cortex", "superior frontal gyrus", "orbitofrontal cortex",
        "medial orbital frontal cortex", "Broca's area", "prefrontal cortex",
        "temporal lobe", "left temporal lobe", "right temporal lobe", 
        "angular gyrus", "entorhinal cortex",
        "parietal lobe", "left parietal lobe", "right parietal lobe", "primary somatosensory cortex",
        "occipital lobe", "right occipital lobe", "primary visual cortex",
        "occipital cortex", "insular cortex", "parietal cortex", "temporal cortex",
        "frontal cortex", "Brodmann (1909) area 4", "temporoparietal junction",
        "middle temporal gyrus", "cingulate cortex", "brain", "brain white matter", "cerebral cortex", "cerebral nuclei"
      ),
      
      # Limbic and Basal Systems
      "limbic and basal systems" = c(
        "anterior cingulate cortex", "anterior cingulate gyrus", "hippocampal formation",
        "hypothalamus", "thalamic complex", "dentate nucleus", "basal ganglion",
        "caudate nucleus", "putamen", "substantia nigra pars compacta",
        "lateral ganglionic eminence", "medial ganglionic eminence",
        "caudal ganglionic eminence", "ganglionic eminence"
      ),
      
      # Brainstem and Cerebellar Structures
      "brainstem and cerebellar structures" = c(
        "pons", "midbrain", "myelencephalon", "telencephalon", "forebrain",
        "cerebellum", "cerebellum vermis lobule", "cerebellar cortex",
        "hemisphere part of cerebellar posterior lobe", "white matter of cerebellum"
      ),
      
      # General Brain and Major Structures
      "general brain and major structures" = c(
        "spinal cord", "neural tube", "cervical spinal cord white matter"
      ),
      
      # Muscular System (Skeletal Muscles)
      "muscular system (skeletal muscles)" = c(
        "rectus abdominis muscle", "gastrocnemius", "muscle of abdomen", "muscle organ",
        "muscle tissue", "pelvic diaphragm muscle", "skeletal muscle tissue", "muscle of pelvic diaphragm"
      ),
      
      # Connective Tissue
      "connective tissue" = c(
        "connective tissue", "tendon of semitendinosus", "vault of skull", "bone spine",
        "rib"
      ),
      
      # Adipose Tissue
      "adipose tissue" = c(
        "adipose tissue", "subcutaneous adipose tissue", "visceral abdominal adipose tissue",
        "perirenal fat", "omental fat pad", "subcutaneous abdominal adipose tissue",
        "abdominal adipose tissue"
      ),
      
      # Endocrine System
      "endocrine system" = c(
        "thyroid gland", "adrenal tissue", "adrenal gland", "islet of Langerhans",
        "endocrine pancreas", "pineal gland"
      ),
      
      # Lymphatic System
      "lymphatic system" = c(
        "lymph node", "mesenteric lymph node", "thoracic lymph node",
        "cervical lymph node", "bronchopulmonary lymph node", "tonsil", "inguinal lymph node"
      ),
      
      # Integumentary System (Skin)
      "integumentary system (skin)" = c(
        "skin of abdomen", "skin of forearm", "skin of scalp", "skin of face", "skin of leg",
        "skin of chest", "skin of back", "skin of hip", "skin of body", "skin of cheek",
        "skin of temple", "skin of shoulder", "skin of external ear", "skin of trunk",
        "skin of prepuce of penis", "skin epidermis", "arm skin", "lower leg skin",
        "hindlimb skin", "zone of skin", "dermis", "skin of nose", "skin of forehead",
        "skin of pes", "axilla"
      ),
      
      # Gastrointestinal Accessory Organs
      "gallbladder" =  "gallbladder",
      
      # Gastrointestinal Accessory Organs
      "pancreas" = c( "pancreas", "exocrine pancreas" ),
      
      # Gastrointestinal Accessory Organs
      "liver" = c( "liver", "caudate lobe of liver", "hepatic cecum" ),
      
      # Spleen
      "spleen" = "spleen",
      
      # Thymus
      "thymus" = "thymus",
      
      # Blood
      "blood" = "blood",
      
      # Bone Marrow
      "bone marrow" = "bone marrow",
      
      # Female Reproductive System
      "female reproductive system" = c(
        "uterus", "myometrium", "fallopian tube", "ampulla of uterine tube",
        "fimbria of uterine tube", "uterine cervix", "endometrium",
        "decidua", "decidua basalis", "placenta", "yolk sac", "isthmus of fallopian tube"
      ),
      "ovary" = "ovary", 
      
      # Male Reproductive System
      "male reproductive system (other)" = c(
        "testis", "gonad"
      ),
      
      # Prostate
      "prostate" = c(
        "prostate gland", "transition zone of prostate", "peripheral zone of prostate", "prostate"
      ),
      
      # Renal System
      "renal system" = c(
        "kidney", "cortex of kidney", "renal medulla", "renal papilla",
        "renal pelvis", "ureter", "bladder organ"
      ),
      
      # Miscellaneous Glands
      "miscellaneous glands" = c(
        "parotid gland", "lacrimal gland", "sublingual gland", "mammary gland",
        "chorionic villus"
      ),
      
      # Eye and Visual-Related Structures
      "sensory-related structures" = c(
        "retina",
        "retinal neural layer",
        "macula lutea",
        "macula lutea proper",
        "sclera",
        "trabecular meshwork",
        "conjunctiva",
        "pigment epithelium of eye",
        "cornea",
        "iris",
        "ciliary body",
        "peripheral region of retina",
        "eye trabecular meshwork",
        "perifoveal part of retina",
        "choroid plexus",
        "lens of camera-type eye",
        "corneo-scleral junction",
        "fovea centralis",
        "eye",
        "inner ear",
        "vestibular system",
        "primary auditory cortex"
      ),
      
      # Digestive Tract Junctions and Connections
      "digestive tract junctions and connections" = c(
        "esophagogastric junction", "duodeno-jejunal junction", "hepatopancreatic ampulla",
        "hepatopancreatic duct", "pyloric antrum"
      ),
      
      # Peritoneal and Abdominal Cavity Structures
      "peritoneal and abdominal cavity structures" = c(
        "peritoneum", "omentum", "retroperitoneum", "mesentery"
      ),
      
      # Breast
      "breast" = c(
        "breast", "upper outer quadrant of breast"
      )
    )
    
    # Add HCA2025
    
    # Respiratory system
    tissue_grouped[["respiratory system"]] <- unique(c(
      tissue_grouped[["respiratory system"]],
      "epithelium of main bronchus",
      "alveolar sac",
      "terminal bronchus",
      "segmental bronchus",
      "lung epithelium",
      "epithelial lining fluid",
      "bronchial mucosa"
    ))
    
    # Trachea
    tissue_grouped[["trachea"]] <- unique(c(
      tissue_grouped[["trachea"]],
      "carina of trachea"
    ))
    
    # Cardiovascular system
    tissue_grouped[["cardiovascular system"]] <- unique(c(
      tissue_grouped[["cardiovascular system"]],
      "right coronary artery",
      "ascending aorta",
      "descending aorta",
      "bulb of aorta",
      "arch of aorta",
      "abdominal aorta",
      "right atrium auricular region",
      "left atrium auricular region",
      "apical region of left ventricle"
    ))
    
    # Vasculature
    tissue_grouped[["vasculature"]] <- unique(c(
      tissue_grouped[["vasculature"]],
      "iliac artery",
      "pulmonary artery",
      "femoral artery",
      "carotid artery segment"
    ))
    
    # Oesophagus
    tissue_grouped[["oesophagus"]] <- unique(c(
      tissue_grouped[["oesophagus"]],
      "esophagus muscularis mucosa"
    ))
    
    # Stomach
    tissue_grouped[["stomach"]] <- unique(c(
      tissue_grouped[["stomach"]],
      "epithelium of stomach",
      "pylorus"
    ))
    
    # Small intestine
    tissue_grouped[["small intestine"]] <- unique(c(
      tissue_grouped[["small intestine"]],
      "ileum lamina propria",
      "jejunum lamina propria",
      "duodenal epithelium"
    ))
    
    # Large intestine
    tissue_grouped[["large intestine"]] <- unique(c(
      tissue_grouped[["large intestine"]],
      "mucosa of sigmoid colon",
      "colorectum",
      "colonic mucosa",
      "cecum mucosa",
      "mucosa of transverse colon",
      "rectosigmoid junction"
    ))
    
    # Digestive system
    tissue_grouped[["digestive system (general)"]] <- unique(c(
      tissue_grouped[["digestive system (general)"]],
      "intestinal mucosa"
    ))
    
    # Nasal, oral and pharyngeal
    tissue_grouped[["nasal, oral, and pharyngeal regions"]] <- unique(c(
      tissue_grouped[["nasal, oral, and pharyngeal regions"]],
      "oropharynx",
      "nasal cavity epithelium",
      "oral cavity",
      "dental pulp"
    ))
    
    # Cerebral lobes / cortex
    tissue_grouped[["cerebral lobes and cortical areas"]] <- unique(c(
      tissue_grouped[["cerebral lobes and cortical areas"]],
      "cortex",
      "neocortex",
      "cortical plate",
      "lower layers of the cortex",
      "upper layers of the cortex",
      "lobe of cerebral hemisphere",
      "Brodmann (1909) area 9",
      "Brodmann (1909) area 10",
      "Brodmann (1909) area 7",
      "Brodmann (1909) area 46",
      "Brodmann (1909) area 17",
      "white matter of temporal lobe",
      "white matter of occipital lobe",
      "white matter of parietal lobe",
      "white matter of frontal lobe",
      "visual cortex"
    ))
    
    # Limbic and basal systems
    tissue_grouped[["limbic and basal systems"]] <- unique(c(
      tissue_grouped[["limbic and basal systems"]],
      "medial globus pallidus",
      "nucleus accumbens",
      "ventral part of telencephalon",
      "dorsal plus ventral thalamus",
      "cerebral subcortex",
      "diencephalon",
      "basal forebrain",
      "arcuate nucleus of hypothalamus",
      "amygdala",
      "striatum",
      "caudate-putamen"
    ))
    
    # Brainstem / cerebellar structures
    tissue_grouped[["brainstem and cerebellar structures"]] <- unique(c(
      tissue_grouped[["brainstem and cerebellar structures"]],
      "metencephalon",
      "hindbrain",
      "dorsal motor nucleus of vagus nerve",
      "cerebellar peduncular complex",
      "midbrain tectum",
      "medulla oblongata",
      "midbrain tegmentum"
    ))
    
    # General brain and major structures
    tissue_grouped[["general brain and major structures"]] <- unique(c(
      tissue_grouped[["general brain and major structures"]],
      "brain meninx",
      "forebrain meninges",
      "meninx",
      "cerebrospinal fluid"
    ))
    
    # Muscular system
    tissue_grouped[["muscular system (skeletal muscles)"]] <- unique(c(
      tissue_grouped[["muscular system (skeletal muscles)"]],
      "intercostal muscle",
      "skeletal muscle organ, vertebrate",
      "diaphragm"
    ))
    
    # Connective tissue
    tissue_grouped[["connective tissue"]] <- unique(c(
      tissue_grouped[["connective tissue"]],
      "calcaneal tendon",
      "tendon of quadriceps femoris"
    ))
    
    # Adipose
    tissue_grouped[["adipose tissue"]] <- unique(c(
      tissue_grouped[["adipose tissue"]],
      "visceral fat"
    ))
    
    # Lymphatic system
    tissue_grouped[["lymphatic system"]] <- unique(c(
      tissue_grouped[["lymphatic system"]],
      "external iliac lymph node",
      "cystic lymph node"
    ))
    
    # Skin
    tissue_grouped[["integumentary system (skin)"]] <- unique(c(
      tissue_grouped[["integumentary system (skin)"]],
      "skin of neck",
      "skin of knee",
      "nose skin",
      "scalp",
      "skin of calf"
    ))
    
    # Liver
    tissue_grouped[["liver"]] <- unique(c(
      tissue_grouped[["liver"]],
      "intrahepatic bile duct",
      "left lobe of liver",
      "right lobe of liver"
    ))
    
    # Gallbladder / biliary system
    tissue_grouped[["gallbladder"]] <- unique(c(
      tissue_grouped[["gallbladder"]],
      "biliary system"
    ))
    
    # Female reproductive system
    tissue_grouped[["female reproductive system"]] <- unique(c(
      tissue_grouped[["female reproductive system"]],
      "fimbria of fallopian tube",
      "trophoblast",
      "adnexa of uterus",
      "ampulla of fallopian tube"
    ))
    
    # Ovary
    tissue_grouped[["ovary"]] <- unique(c(
      tissue_grouped[["ovary"]],
      "left ovary",
      "right ovary"
    ))
    
    # Renal / urinary system
    tissue_grouped[["renal system"]] <- unique(c(
      tissue_grouped[["renal system"]],
      "urinary bladder"
    ))
    
    # Miscellaneous glands
    tissue_grouped[["miscellaneous glands"]] <- unique(c(
      tissue_grouped[["miscellaneous glands"]],
      "labial gland",
      "submandibular gland",
      "salivary gland epithelium"
    ))
    
    # Sensory / ocular structures
    tissue_grouped[["sensory-related structures"]] <- unique(c(
      tissue_grouped[["sensory-related structures"]],
      "optic choroid",
      "corneal epithelium",
      "corneal endothelium",
      "chorioretinal region",
      "substantia propria of cornea",
      "uvea",
      "choroid plexus epithelium",
      "cranial nerve II",
      "ocular surface region",
      "pigmented layer of retina",
      "tunica fibrosa of eyeball",
      "eyeball of camera-type eye",
      "optic disc",
      "central retina",
      "dorsal root ganglion"
    ))
    
    # Digestive tract junctions / connections
    tissue_grouped[["digestive tract junctions and connections"]] <- unique(c(
      tissue_grouped[["digestive tract junctions and connections"]],
      "common bile duct"
    ))
    
    # Peritoneal / abdominal cavity
    tissue_grouped[["peritoneal and abdominal cavity structures"]] <- unique(c(
      tissue_grouped[["peritoneal and abdominal cavity structures"]],
      "paracolic gutter",
      "parietal peritoneum",
      "ascitic fluid"
    ))
    
    # Breast
    tissue_grouped[["breast"]] <- unique(c(
      tissue_grouped[["breast"]],
      "mammary gland connective tissue",
      "parenchyma of mammary gland",
      "milk"
    ))
    
    # Reshape
    tissue_grouped = 
      tissue_grouped |> 
      enframe(name ="tissue_groups") |> 
      distinct() |> 
      unnest(value) |> 
      rename(tissue = value) 
    
    ethnicity_grouped <- tribble(
      ~self_reported_ethnicity, ~ethnicity_groups,
      "unknown", "Other/Unknown",
      "European", "European",
      "Korean", "East Asian",
      "Asian", "East Asian",
      "Japanese", "Japanese",
      "African American", "African",
      "Hispanic or Latin American", "Hispanic/Latin American",
      "Hispanic or Latin", "Hispanic/Latin American",
      "Singaporean Chinese", "East Asian",
      "Han Chinese", "East Asian",
      "Singaporean Indian", "South Asian",
      "Singaporean Malay", "Other/Unknown",
      "British", "European",
      "African", "African",
      "South Asian", "South Asian",
      "European American", "European",
      "East Asian", "East Asian",
      "American", "Other/Unknown",
      "African American or Afro-Caribbean", "African",
      "Oceanian", "Native American & Pacific Islander",
      "Jewish Israeli", "Middle Eastern & North African",
      "Chinese", "East Asian",
      "South East Asian", "Other/Unknown",
      "Greater Middle Eastern  (Middle Eastern or North African or Persian)", "Middle Eastern & North African",
      "Native American", "Native American & Pacific Islander",
      "Pacific Islander", "Native American & Pacific Islander",
      "Finnish", "European",
      "Bangladeshi", "South Asian",
      "Native American,Hispanic or Latin American", "Hispanic/Latin American",
      "Irish", "European",
      "Iraqi", "Middle Eastern & North African",
      "European,Asian", "European",
      "White", "European",
      "Hispanic", "Hispanic/Latin American"
    )
    
    ethnicity_grouped <- bind_rows(
      ethnicity_grouped,
      tribble(
        ~self_reported_ethnicity, ~ethnicity_groups,
        "Australian", "Other/Unknown",
        "Canadian", "Other/Unknown",
        "Thai", "East Asian",
        "European American || Asian", "European",
        "Ethiopian", "African",
        "Arab", "Middle Eastern & North African",
        "Southeast Asian", "East Asian",
        "Indian", "South Asian",
        "Brazilian", "Hispanic/Latin American",
        "Hispanic or Latin || Native American", "Hispanic/Latin American",
        "European American || Hispanic or Latin", "Hispanic/Latin American",
        "African American || European American || Native American", "Other/Unknown",
        "German", "European",
        "Middle Eastern", "Middle Eastern & North African"
      )
    ) %>%
      distinct(self_reported_ethnicity, .keep_all = TRUE)
    
    assay_data_grouped <- tribble(
      ~assay, ~assay_groups,
      "10x 3' v2", "10x Genomics 3",
      "10x 3' v3", "10x Genomics 3",
      "10x 5' v2", "10x Genomics 5",
      "10x 5' v1", "10x Genomics 5",
      "10X Genomics Chromium", "10x Genomics",
      "MARS-seq", "Plate based Technologies",
      "10x 3' transcription profiling", "10x Genomics 3",
      "10x 5' transcription profiling", "10x Genomics 5",
      "Smart-seq2", "Smart seq",
      "microwell-seq", "Microwell Technologies",
      "TruDrop", "TruDrop",
      "Drop-seq", "Drop based Technologies",
      "Seq-Well S3", "Microwell Technologies",
      "GEXSCOPE technology", "Other Technologies",
      "Seq-Well", "Microwell Technologies",
      "sci-RNA-seq", "Other Technologies",
      "10x 3' v1", "10x Genomics 3",
      "BD Rhapsody Whole Transcriptome Analysis", "Other Technologies",
      "BD Rhapsody Targeted mRNA", "Other Technologies",
      "CEL-seq2", "Plate based Technologies",
      "SPLiT-seq", "Other Technologies",
      "STRT-seq", "Plate based Technologies",
      "inDrop", "Drop based Technologies",
      "Smart-seq v4", "Smart seq",
      "ScaleBio single cell RNA sequencing", "Other Technologies"
    )
    
    assay_data_grouped <- bind_rows(
      assay_data_grouped,
      tribble(
        ~assay, ~assay_groups,
        "mCT-seq", "Other Technologies",
        "10x multiome", "10x Genomics 3",
        "SORT-seq", "Plate based Technologies",
        "Parse Evercode Whole Transcriptome v2", "Other Technologies",
        "particle-templated instant partition sequencing", "Drop based Technologies",
        "Quartz-seq", "Plate based Technologies",
        "CEL-seq", "Plate based Technologies",
        "Asteria scRNA-seq kit", "Other Technologies",
        "HIVE CLX Single-Cell RNAseq Solution", "Microwell Technologies",
        "10x gene expression flex", "10x Genomics 3",
        "modified STRT-seq", "Plate based Technologies",
        "Smart-seq", "Smart seq",
        "sci-RNA-seq3", "Other Technologies"
      )
    ) %>%
      distinct(assay, .keep_all = TRUE)
    
    disease_data_grouped <- tribble(
      ~disease, ~disease_groups,
      
      # Normal control
      "normal", "Normal",
      
      # Isolated Diseases
      "COVID-19", "COVID-19 related",
      "post-COVID-19 disorder", "COVID-19 related",
      "long COVID-19", "COVID-19 related",
      "glioblastoma", "Glioblastoma",
      "lung adenocarcinoma", "Lung Adenocarcinoma",
      "systemic lupus erythematosus", "Systemic Lupus Erythematosus",
      
      # Infectious and Immune-related Diseases (other than COVID-19)
      "Crohn disease", "Infectious and Immune-related Diseases",
      "Crohn ileitis", "Infectious and Immune-related Diseases",
      "pneumonia", "Infectious and Immune-related Diseases",
      "common variable immunodeficiency", "Infectious and Immune-related Diseases",
      "toxoplasmosis", "Infectious and Immune-related Diseases",
      "Plasmodium malariae malaria", "Infectious and Immune-related Diseases",
      "type 1 diabetes mellitus", "Infectious and Immune-related Diseases",
      "influenza", "Infectious and Immune-related Diseases",
      "chronic rhinitis", "Infectious and Immune-related Diseases",
      "periodontitis", "Infectious and Immune-related Diseases",
      "localized scleroderma", "Infectious and Immune-related Diseases",
      "lymphangioleiomyomatosis", "Infectious and Immune-related Diseases",
      "listeriosis", "Infectious and Immune-related Diseases",
      
      # Cancer (other than isolated cancers)
      "squamous cell lung carcinoma", "Cancer",
      "small cell lung carcinoma", "Cancer",
      "non-small cell lung carcinoma", "Cancer",
      "breast carcinoma", "Cancer",
      "breast cancer", "Cancer",
      "luminal B breast carcinoma", "Cancer",
      "luminal A breast carcinoma", "Cancer",
      "triple-negative breast carcinoma", "Cancer",
      "gastric cancer", "Cancer",
      "colorectal cancer", "Cancer",
      "colon sessile serrated adenoma/polyp", "Cancer",
      "follicular lymphoma", "Cancer",
      "B-cell acute lymphoblastic leukemia", "Cancer",
      "B-cell non-Hodgkin lymphoma", "Cancer",
      "acute myeloid leukemia", "Cancer",
      "acute promyelocytic leukemia", "Cancer",
      "plasma cell myeloma", "Cancer",
      "clear cell renal carcinoma", "Cancer",
      "nonpapillary renal cell carcinoma", "Cancer",
      "basal cell carcinoma", "Cancer",
      "colorectal neoplasm", "Cancer",
      "adenocarcinoma", "Cancer",
      "chromophobe renal cell carcinoma", "Cancer",
      "neuroendocrine carcinoma", "Cancer",
      "lung large cell carcinoma", "Cancer",
      "tongue cancer", "Cancer",
      "Wilms tumor", "Cancer",
      "pleomorphic carcinoma", "Cancer",
      "blastoma", "Cancer",
      
      # Neurodegenerative and Neurological Disorders
      "dementia", "Neurodegenerative and Neurological Disorders",
      "Alzheimer disease", "Neurodegenerative and Neurological Disorders",
      "Parkinson disease", "Neurodegenerative and Neurological Disorders",
      "amyotrophic lateral sclerosis", "Neurodegenerative and Neurological Disorders",
      "multiple sclerosis", "Neurodegenerative and Neurological Disorders",
      "Down syndrome", "Neurodegenerative and Neurological Disorders",
      "trisomy 18", "Neurodegenerative and Neurological Disorders",
      "frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
      "temporal lobe epilepsy", "Neurodegenerative and Neurological Disorders",
      "Lewy body dementia", "Neurodegenerative and Neurological Disorders",
      "amyotrophic lateral sclerosis 26 with or without frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
      
      # Respiratory Conditions
      "pulmonary fibrosis", "Respiratory Conditions",
      "respiratory system disorder", "Respiratory Conditions",
      "chronic obstructive pulmonary disease", "Respiratory Conditions",
      "cystic fibrosis", "Respiratory Conditions",
      "interstitial lung disease", "Respiratory Conditions",
      "hypersensitivity pneumonitis", "Respiratory Conditions",
      "non-specific interstitial pneumonia", "Respiratory Conditions",
      "aspiration pneumonia", "Respiratory Conditions",
      "pulmonary emphysema", "Respiratory Conditions",
      "pulmonary sarcoidosis", "Respiratory Conditions",
      
      # Cardiovascular Diseases
      "myocardial infarction", "Cardiovascular Diseases",
      "acute myocardial infarction", "Cardiovascular Diseases",
      "dilated cardiomyopathy", "Cardiovascular Diseases",
      "heart failure", "Cardiovascular Diseases",
      "arrhythmogenic right ventricular cardiomyopathy", "Cardiovascular Diseases",
      "congenital heart disease", "Cardiovascular Diseases",
      "non-compaction cardiomyopathy", "Cardiovascular Diseases",
      "cardiomyopathy", "Cardiovascular Diseases",
      "heart disorder", "Cardiovascular Diseases",
      
      # Metabolic and Other Disorders
      "type 2 diabetes mellitus", "Metabolic and Other Disorders",
      "chronic kidney disease", "Metabolic and Other Disorders",
      "digestive system disorder", "Metabolic and Other Disorders",
      "primary sclerosing cholangitis", "Metabolic and Other Disorders",
      "gastritis", "Metabolic and Other Disorders",
      "acute kidney failure", "Metabolic and Other Disorders",
      "tubular adenoma", "Metabolic and Other Disorders",
      "benign prostatic hyperplasia", "Metabolic and Other Disorders",
      "opiate dependence", "Metabolic and Other Disorders",
      "gingivitis", "Metabolic and Other Disorders",
      "hyperplastic polyp", "Metabolic and Other Disorders",
      "clonal hematopoiesis", "Metabolic and Other Disorders",
      "epilepsy", "Metabolic and Other Disorders",
      "age related macular degeneration 7", "Metabolic and Other Disorders",
      "kidney benign neoplasm", "Metabolic and Other Disorders",
      "malignant pancreatic neoplasm", "Metabolic and Other Disorders",
      "cataract", "Metabolic and Other Disorders",
      "macular degeneration", "Metabolic and Other Disorders",
      "hydrosalpinx", "Metabolic and Other Disorders",
      "tubulovillous adenoma", "Metabolic and Other Disorders",
      "gastric intestinal metaplasia", "Metabolic and Other Disorders",
      "Barrett esophagus", "Metabolic and Other Disorders",
      
      # Other Diseases
      "injury", "Other Diseases",
      "anencephaly", "Other Diseases",
      "primary biliary cholangitis", "Other Diseases",
      "keloid", "Other Diseases",
      "kidney oncocytoma", "Other Diseases",
      "respiratory failure", "Other Diseases",
      "pilocytic astrocytoma", "Other Diseases"
    )
    
    new_disease_data_grouped <- tribble(
      ~disease, ~disease_groups,
      
      "dementia || Parkinson disease", "Neurodegenerative and Neurological Disorders",
      "malignant ovarian serous tumor", "Cancer",
      "leukoencephalopathy, diffuse hereditary, with spheroids 1", "Neurodegenerative and Neurological Disorders",
      "dementia || progressive supranuclear palsy", "Neurodegenerative and Neurological Disorders",
      "Sjogren syndrome", "Infectious and Immune-related Diseases",
      "vascular dementia || Lewy body dementia", "Neurodegenerative and Neurological Disorders",
      "schizophrenia || type 1 diabetes mellitus", "Other Diseases",
      "breast mucinous carcinoma", "Cancer",
      "dementia || Alzheimer disease", "Neurodegenerative and Neurological Disorders",
      "dementia || amyotrophic lateral sclerosis", "Neurodegenerative and Neurological Disorders",
      "tauopathy", "Neurodegenerative and Neurological Disorders",
      "oral cavity squamous cell carcinoma", "Cancer",
      "melanoma", "Cancer",
      "dementia || post-traumatic stress disorder", "Neurodegenerative and Neurological Disorders",
      "dementia || schizophrenia || head injury", "Neurodegenerative and Neurological Disorders",
      "bipolar I disorder || brain neoplasm", "Other Diseases",
      "dementia || normal pressure hydrocephalus", "Neurodegenerative and Neurological Disorders",
      "dementia || Parkinson disease || Lewy body dementia || head injury", "Neurodegenerative and Neurological Disorders",
      "atrial fibrillation || mitral valve insufficiency", "Cardiovascular Diseases",
      "open-angle glaucoma", "Metabolic and Other Disorders",
      "prediabetes syndrome", "Metabolic and Other Disorders",
      "dementia || vascular dementia || Alzheimer disease || schizophrenia", "Neurodegenerative and Neurological Disorders",
      "hematologic disorder", "Other Diseases",
      "isolated focal cortical dysplasia type II", "Neurodegenerative and Neurological Disorders",
      "type 1 diabetes mellitus || schizoaffective depressive disorder", "Other Diseases",
      "bipolar I disorder || bulimia nervosa", "Other Diseases",
      "dementia || Alzheimer disease || normal pressure hydrocephalus", "Neurodegenerative and Neurological Disorders",
      "bipolar I disorder || attention deficit-hyperactivity disorder", "Other Diseases",
      "schizophrenia || type 1 diabetes mellitus || tardive dyskinesia", "Other Diseases",
      "bipolar I disorder || atherosclerosis", "Other Diseases",
      "estrogen-receptor positive breast cancer", "Cancer",
      "vascular dementia || progressive supranuclear palsy", "Neurodegenerative and Neurological Disorders",
      "dementia || vascular dementia || diabetes mellitus", "Neurodegenerative and Neurological Disorders",
      "schizophrenia || type 2 diabetes mellitus || tardive dyskinesia", "Other Diseases",
      "major depressive disorder || schizophrenia || atherosclerosis", "Other Diseases",
      "vascular dementia || diabetes mellitus || Parkinson disease", "Neurodegenerative and Neurological Disorders",
      "dementia || major depressive disorder || Alzheimer disease || Lewy body dementia", "Neurodegenerative and Neurological Disorders",
      "bipolar I disorder || type 2 diabetes mellitus", "Other Diseases",
      "brain cancer", "Cancer",
      "bipolar I disorder", "Other Diseases",
      "dementia || vascular dementia", "Neurodegenerative and Neurological Disorders",
      "myelodysplastic/myeloproliferative disease", "Cancer",
      "bipolar I disorder || anorexia nervosa || obsessive-compulsive disorder", "Other Diseases",
      "invasive tubular breast carcinoma || invasive lobular breast carcinoma", "Cancer",
      "schizophrenia || brain neoplasm", "Other Diseases",
      "HER2 positive breast carcinoma", "Cancer",
      "bipolar disorder", "Other Diseases",
      "type 1 diabetes mellitus || atherosclerosis", "Other Diseases",
      "type 2 diabetes mellitus || atherosclerosis", "Other Diseases",
      "cytomegalovirus infection", "Infectious and Immune-related Diseases",
      "dementia || Alzheimer disease || Lewy body dementia", "Neurodegenerative and Neurological Disorders",
      "myelodysplastic syndrome", "Cancer",
      "cognitive disorder", "Neurodegenerative and Neurological Disorders",
      "dementia || schizophrenia || Parkinson disease", "Neurodegenerative and Neurological Disorders",
      "dementia || frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
      "dementia || Alzheimer disease || brain neoplasm", "Neurodegenerative and Neurological Disorders",
      "dementia || diabetes mellitus || Parkinson disease", "Neurodegenerative and Neurological Disorders",
      "HIV infectious disease || leishmaniasis", "Infectious and Immune-related Diseases",
      "head injury", "Other Diseases",
      "vascular dementia || diabetes mellitus", "Neurodegenerative and Neurological Disorders",
      "juvenile dermatomyositis", "Infectious and Immune-related Diseases",
      "bipolar I disorder || obsessive-compulsive disorder", "Other Diseases",
      "atherosclerosis || tardive dyskinesia || schizoaffective bipolar disorder", "Other Diseases",
      "squamous cell carcinoma", "Cancer",
      "enamel caries", "Other Diseases",
      "brain neoplasm", "Cancer",
      "dementia || vascular dementia || Alzheimer disease || diabetes mellitus", "Neurodegenerative and Neurological Disorders",
      "bipolar I disorder || post-traumatic stress disorder", "Other Diseases",
      "bipolar II disorder || type 2 diabetes mellitus", "Other Diseases",
      "dementia || Lewy body dementia", "Neurodegenerative and Neurological Disorders",
      "atherosclerosis", "Cardiovascular Diseases",
      "dementia || vascular dementia || tauopathy", "Neurodegenerative and Neurological Disorders",
      "dementia || multiple sclerosis", "Neurodegenerative and Neurological Disorders",
      "bipolar I disorder || type 1 diabetes mellitus", "Other Diseases",
      "bipolar II disorder || type 1 diabetes mellitus", "Other Diseases",
      "heart valve disorder", "Cardiovascular Diseases",
      "bipolar I disorder || anorexia nervosa", "Other Diseases",
      "lung cancer", "Cancer",
      "myocarditis", "Cardiovascular Diseases",
      "lymphoid interstitial pneumonia", "Respiratory Conditions",
      "dementia || vascular dementia || head injury", "Neurodegenerative and Neurological Disorders",
      "dementia || brain neoplasm", "Neurodegenerative and Neurological Disorders",
      "normal pressure hydrocephalus", "Neurodegenerative and Neurological Disorders",
      "schizophrenia || type 1 diabetes mellitus || atherosclerosis", "Other Diseases",
      "ulcerative colitis", "Infectious and Immune-related Diseases",
      "ovarian cancer", "Cancer",
      "major depressive disorder || tauopathy", "Neurodegenerative and Neurological Disorders",
      "vascular dementia", "Neurodegenerative and Neurological Disorders",
      "oropharynx squamous cell carcinoma", "Cancer",
      "schizoaffective depressive disorder", "Other Diseases",
      "dementia || Alzheimer disease || Parkinson disease || Lewy body dementia", "Neurodegenerative and Neurological Disorders",
      "epilepsy || schizophrenia || tardive dyskinesia", "Other Diseases",
      "metastatic melanoma", "Cancer",
      "epilepsy || schizophrenia", "Other Diseases",
      "lymphoma", "Cancer",
      "dementia || Alzheimer disease || diabetes mellitus", "Neurodegenerative and Neurological Disorders",
      "dementia || Alzheimer disease || Parkinson disease", "Neurodegenerative and Neurological Disorders",
      "anthracosis", "Respiratory Conditions",
      "coronary artery disorder", "Cardiovascular Diseases",
      "HIV infectious disease || visceral leishmaniasis", "Infectious and Immune-related Diseases",
      "hypertrophic cardiomyopathy", "Cardiovascular Diseases",
      "primary cutaneous diffuse large B-cell lymphoma, Leg type", "Cancer",
      "dementia || Parkinson disease || progressive supranuclear palsy", "Neurodegenerative and Neurological Disorders",
      "diabetes mellitus", "Metabolic and Other Disorders",
      "neuroblastoma", "Cancer",
      "bipolar I disorder || bulimia nervosa || obsessive-compulsive disorder", "Other Diseases",
      "dementia || tauopathy", "Neurodegenerative and Neurological Disorders",
      "schizophrenia || atherosclerosis", "Other Diseases",
      "intestinal failure–associated liver disease", "Metabolic and Other Disorders",
      "diffuse large B-cell lymphoma", "Cancer",
      "celiac disease", "Infectious and Immune-related Diseases",
      "dementia || Alzheimer disease || amyotrophic lateral sclerosis", "Neurodegenerative and Neurological Disorders",
      "dementia || Alzheimer disease || schizophrenia || Lewy body dementia", "Neurodegenerative and Neurological Disorders",
      "schizophrenia || type 2 diabetes mellitus", "Other Diseases",
      "myocardial ischemia", "Cardiovascular Diseases",
      "colon adenocarcinoma", "Cancer",
      "prostatic acinar adenocarcinoma", "Cancer",
      "inflammatory bowel disease", "Infectious and Immune-related Diseases",
      "dementia || head injury", "Neurodegenerative and Neurological Disorders",
      "bipolar II disorder", "Other Diseases",
      "pulpitis", "Infectious and Immune-related Diseases",
      "dementia || vascular dementia || Alzheimer disease || Parkinson disease", "Neurodegenerative and Neurological Disorders",
      "bipolar disorder || attention deficit-hyperactivity disorder", "Other Diseases",
      "schizophrenia", "Other Diseases",
      "dementia || vascular dementia || schizophrenia", "Neurodegenerative and Neurological Disorders",
      "bronchopulmonary dysplasia", "Respiratory Conditions",
      "dementia || Alzheimer disease || schizophrenia", "Neurodegenerative and Neurological Disorders",
      "progressive supranuclear palsy", "Neurodegenerative and Neurological Disorders",
      "colorectal carcinoma || metastatic malignant neoplasm", "Cancer",
      "rheumatoid arthritis", "Infectious and Immune-related Diseases",
      "liver cancer", "Cancer",
      "diabetes mellitus || schizophrenia", "Other Diseases",
      "inflammatory disease", "Infectious and Immune-related Diseases",
      "dementia || vascular dementia || Alzheimer disease", "Neurodegenerative and Neurological Disorders",
      "dementia || schizophrenia", "Neurodegenerative and Neurological Disorders",
      "dementia || Alzheimer disease || frontotemporal dementia", "Neurodegenerative and Neurological Disorders",
      "diabetes mellitus || head injury", "Other Diseases",
      "dementia || Alzheimer disease || diabetes mellitus || Parkinson disease || head injury",
      "Neurodegenerative and Neurological Disorders",
      "obstructive nephropathy", "Metabolic and Other Disorders",
      "schizophrenia || tardive dyskinesia", "Other Diseases",
      "Pick disease", "Neurodegenerative and Neurological Disorders",
      "uveal melanoma", "Cancer",
      "HIV infectious disease", "Infectious and Immune-related Diseases",
      "bipolar disorder || schizoaffective bipolar disorder", "Other Diseases",
      "bipolar disorder || type 2 diabetes mellitus || schizoaffective bipolar disorder",
      "Other Diseases"
    )
    
    disease_data_grouped <- bind_rows(
      disease_data_grouped,
      tribble(
        ~disease, ~disease_groups,
        "invasive lobular breast carcinoma", "Cancer",
        "invasive ductal breast carcinoma", "Cancer"
      )
    )
    
    # Add only genuinely new disease strings.
    # Existing mappings always win and therefore cannot be overwritten.
    disease_data_grouped <- bind_rows(
      disease_data_grouped,
      new_disease_data_grouped %>%
        anti_join(disease_data_grouped, by = "disease")
    )
    
    temp_path = tempdir()
    file.copy("~/git_control/cellNexus_development/disease_data_grouped_further.csv", temp_path)
    #system(glue("/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone copy box:/minh_immune_map_disease/disease_data_grouped_further.csv {temp_path}/"))
    
    disease_data_grouped = 
      disease_data_grouped |> 
      left_join(
        read_csv(glue("{temp_path}/disease_data_grouped_further.csv")) |> 
          rename(disease_groups_further = disease_groups)
      ) |> 
      mutate(disease_groups = if_else(!disease_groups_further |> is.na(), disease_groups_further, disease_groups)) |> 
      select(disease,  disease_groups)
    
    tbl |> 
      
      # TISSUE
      select(-any_of("tissue_groups")) |> 
      left_join(tissue_grouped, copy=TRUE) |> 
      
      # TECH
      left_join(assay_data_grouped, copy=TRUE) |> 
      mutate(assay_groups = if_else(is.na(assay), "Other Technologies", assay_groups)) |> 
    
      # DISEASE
      left_join(disease_data_grouped, copy=TRUE) |> 
      
      # make disease tissue specific, omit Normal
      mutate(disease_groups = paste(disease_groups, tissue_groups, sep = "_")) |> 
      mutate(disease_groups = if_else(disease_groups |> str_detect("Normal_.+"), "Normal", disease_groups))  |> 
      
      
      # TEMPORARY. de-group pancreas and liver
      mutate(tissue_groups = case_when(
        
        tissue %in% c("gallbladder") ~ "gallbladder",
        tissue %in% c("pancreas", "exocrine pancreas") ~ "pancreas",
        tissue %in% c("liver", "caudate lobe of liver", "hepatic cecum" ) ~ "liver",
        TRUE ~ tissue_groups
      )) |> 
      
      # SEX edit
      mutate(sex = if_else(sex |> is.na(), "unknown", sex)) |> 
      
      # Age
      mutate(age_years = age_days / 365) |> 
      mutate(age_bin = dplyr::case_when(
        age_years < 3 ~ "Infancy",
        age_years < 12 ~ "Childhood",
        age_years < 20 ~ "Adolescence",
        age_years < 40 ~ "Young Adulthood",
        age_years < 50 ~ "Middle Age",
        age_years < 60 ~ "Senior_50",
        age_years < 70 ~ "Senior_60",
        age_years >= 70 ~ "Senior_70",
        TRUE ~ NA_character_
      )) |> 
      mutate(age_decade = ceiling(age_years/10) |> as.integer() |> as.character()) |> 
      
      # left_join(age_bin_table, copy=TRUE) |> 
      
      # ETHNICITY
      left_join(ethnicity_grouped, copy=TRUE) |> 
      mutate(ethnicity_groups = if_else(is.na(self_reported_ethnicity), "Other/Unknown", ethnicity_groups)) |> 
      
      dplyr::select(
        cell_id, sample_id, donor_id, dataset_id, title, collection_id, age_days, age_bin, age_decade, sex,
        ethnicity_groups, tissue_groups, tissue, assay_groups, cell_type_unified_ensemble,
        cell_type, disease_groups, is_immune
      ) |>
      
      # Same as scale(age_days, center = 50 * 365), computed in DuckDB, for memory efficiency
      mutate(age_days_centred = age_days - 50 * 365) |>
      mutate(
        age_days_scaled = age_days_centred / sqrt(
          sum(age_days_centred^2, na.rm = TRUE) /
            (sum(as.integer(!is.na(age_days)), na.rm = TRUE) - 1)
        )
      ) |>
      select(-age_days_centred)
      
    
  }     
  
  create_input_cell_counts = function(
    caq_celltype_level_map,
    ethnicity_imputed, count_by = c("sample_id", "donor_id"), 
    result_directory
  ){
    
    count_by <- match.arg(count_by)
    
    metadata_tbl <- get_metadata(
      cloud_metadata = get_metadata_url(c("hca_2024", "hca_2025"))
    )
    
    con <- dbplyr::remote_con(metadata_tbl)
    
    read_census <- function(version) {
      census <- cellxgene.census::open_soma(census_version = version)
      on.exit(census$close(), add = TRUE)
      
      list(
        obs = census$
          get("census_data")$
          get("homo_sapiens")$
          get("obs")$
          read(
            column_names = c(
              "observation_joinid",
              "dataset_id",
              "sex",
              "tissue",
              "assay",
              "disease",
              "self_reported_ethnicity",
              "cell_type",
              "is_primary_data"
            ),
            value_filter = "is_primary_data == 'TRUE'"
          )$
          concat(),
        
        datasets = census$
          get("census_info")$
          get("datasets")$
          read(
            column_names = c(
              "dataset_id",
              "dataset_title",
              "collection_id"
            )
          )$
          concat()
      )
    }
    
    census_2024 <- read_census("2024-07-01")
    census_2025 <- read_census("2025-11-08")
    
    obs_2024_name      <- paste0("census_obs_2024_", count_by)
    datasets_2024_name <- paste0("census_datasets_2024_", count_by)
    obs_2025_name      <- paste0("census_obs_2025_", count_by)
    datasets_2025_name <- paste0("census_datasets_2025_", count_by)
    
    duckdb::duckdb_register_arrow(
      con, obs_2024_name, census_2024$obs
    )
    duckdb::duckdb_register_arrow(
      con, datasets_2024_name, census_2024$datasets
    )
    duckdb::duckdb_register_arrow(
      con, obs_2025_name, census_2025$obs
    )
    duckdb::duckdb_register_arrow(
      con, datasets_2025_name, census_2025$datasets
    )
    
    # All datasets already included in hca_2024.
    old_dataset_ids <- dplyr::tbl(con, datasets_2024_name) |>
      dplyr::select(dataset_id) |>
      dplyr::distinct()
    
    # Only datasets genuinely added in hca_2025.
    new_dataset_ids <- dplyr::tbl(con, datasets_2025_name) |>
      dplyr::select(dataset_id) |>
      dplyr::distinct() |>
      dplyr::anti_join(old_dataset_ids, by = "dataset_id")
    
    # The two sides of this union are now disjoint by dataset_id.
    census_obs_lookup <- dplyr::union_all(
      dplyr::tbl(con, obs_2024_name) |>
        dplyr::select(-is_primary_data) |>
        dplyr::distinct(
          dataset_id,
          observation_joinid,
          .keep_all = TRUE
        ),
      
      dplyr::tbl(con, obs_2025_name) |>
        dplyr::semi_join(new_dataset_ids, by = "dataset_id") |>
        dplyr::select(-is_primary_data) |>
        dplyr::distinct(
          dataset_id,
          observation_joinid,
          .keep_all = TRUE
        )
    )
    
    census_dataset_lookup <- dplyr::union_all(
      dplyr::tbl(con, datasets_2024_name) |>
        dplyr::distinct(dataset_id, .keep_all = TRUE),
      
      dplyr::tbl(con, datasets_2025_name) |>
        dplyr::semi_join(new_dataset_ids, by = "dataset_id") |>
        dplyr::distinct(dataset_id, .keep_all = TRUE)
    ) |>
      dplyr::rename(title = dataset_title)
    
    my_tbl <- metadata_tbl |>
      dplyr::left_join(
        census_obs_lookup,
        by = c("dataset_id", "observation_joinid")
      ) |>
      dplyr::left_join(
        census_dataset_lookup,
        by = "dataset_id"
      ) |>
      keep_quality_cells()
    
    my_tbl <- my_tbl |>
      dplyr::filter(!sample_id %in% c("e6b21291ca8ed20cb62541c2a1c7b2ac", "20b3bf1fa3f5f28eb62130efbe464f94"))
    
    cols_to_keep <- colnames(my_tbl)
    
    prostate_atlas <- get_metadata(local_metadata = "/vast/projects/cellxgene_curated/prostate_atlas/prostate_metadata.v0.3.3.parquet",
                                   cloud_metadata = NULL)
    
    # Combine prostate
    my_tbl <- my_tbl |>
      union_all(prostate_atlas |>
                  filter(is_primary_data == TRUE), copy = TRUE) |>
      select(all_of(cols_to_keep))
    ########################
    
    my_tbl = 
      my_tbl |>
      
      # Add of info to the metadata
      edit_covariates()  |>
      
      # TISSUE
      filter(!tissue_groups %in% c(
        "muscular system (skeletal muscles)",
        "ovary",
        "vasculature",
        "digestive tract junctions and connections",
        "peritoneal and abdominal cavity structures",
        "connective tissue",
        "miscellaneous glands"
      ), 
      !tissue_groups |> is.na()
      ) |> 
      
      # Filter embrios
      filter(age_days > 365) |> 
      
      # NON immune cells
      mutate(cell_type_unified_ensemble = if_else(is_immune, cell_type_unified_ensemble, "non_immune")) |> 
      
      # IMMUNE CELLS
      # filter(is_immune) |> 
      filter(cell_type_unified_ensemble %in% c("non_immune", "cd8 naive", "cd16 mono", "cd4 tcm", "cd4 th17 em", "granulocyte", "cd4 th1/th17 em", "treg", "b memory", "b naive", "nk", "plasma", "cd4 th2 em", "mast", "cd4 th1 em", "cd8 tem", "mait", "tgd", "cdc", "cd4 fh em", "cd4 naive", "nkt", "macrophage", "cd8 tcm", "cd14 mono", "pdc", "ilc")) 
      

    
    # Drop samples that are enriched for a immune cell type
    samples_to_drop =
      flag_enriched_samples_within_immune_compartment(my_tbl, sample_column = count_by) |>
      filter(flag) |>
      pull(all_of(count_by))
    
    my_tbl = my_tbl |> filter(!.data[[count_by]] %in% !!samples_to_drop)

    
    # # Attach plasma subtype
    # left_join(
    #   sce_scored_delayed_plasma |> 
    #     readRDS() |> 
    #     as_tibble()
    # ) |> 
    #         
    # # calculate median per dataset
    # mutate(median_plasma_score = median(TotalScore, na.rm = TRUE), .by = dataset_id) |> 
    # mutate(cell_type_unified_ensemble = 
    #          case_when(
    #            cell_type_unified_ensemble == "plasma" & !TotalScore |> is.na() & TotalScore > median_plasma_score ~ "plasma_long_lived",
    #            cell_type_unified_ensemble == "plasma" & !TotalScore |> is.na() & TotalScore <= median_plasma_score ~ "plasma_short_lived",
    #            TRUE ~ cell_type_unified_ensemble
    #          ))
    
    
    if(count_by=="donor_id")
      my_tbl =  
      my_tbl |> 
      
      # I add dataset_id because there are duplicated donor_id
      # I add technology because there are 71 donors with multiple technologies. It is acceptable in the grand scheme of things of 4K+ donors
      # Same for disease
      mutate(sample_id = paste0(donor_id, dataset_id, tissue_groups, age_days, assay_groups, disease_groups)) 
    
    my_tbl |> 
      dplyr::count(
        sample_id, donor_id, dataset_id, title, collection_id, age_days, age_bin, age_days_scaled, age_decade,
        sex, ethnicity_groups, tissue_groups, tissue, assay_groups, cell_type_unified_ensemble, is_immune,
        disease_groups) |> 
      mutate(n = as.integer(n)) |> 
      as_tibble() |> 
      
      # Add hierarchy of cell types L1, L2, L3
      left_join(
        caq_celltype_level_map |> 
          
          # Plasma hierarchy
          mutate(L4=L3) |>
          bind_rows(tibble(
            cell_type_unified_harmonised = c("plasma_long_lived", "plasma_short_lived"), 
            L0 = "b",   L1 ="plasma",       L2   ="plasma",     L3  = "plasma",    
            L4 =c("plasma_long_lived", "plasma_short_lived")
          )) |> 
          mutate(L4 = case_when(
            L4 |> is.na() ~ NA,
            cell_type_unified_harmonised %in% c("plasma_long_lived", "plasma_short_lived") ~ L4,
            TRUE ~ "non_plasma"
          )),
        
        by = join_by(cell_type_unified_ensemble == cell_type_unified_harmonised)
      ) |> 
      mutate(
        across(
          matches("^L[0-9]"),        # select columns whose names start L0, L1, … L9
          ~ if_else(is_immune, as.character(.), "non_immune")
        )
      ) |> 
      
      # Add imputed ethnicities, and assign original if not present (cmposition and DE might have different samples because of filtering)
      left_join(
        ethnicity_imputed,
        by = join_by(sample_id, ethnicity_groups)
      ) |> 
      mutate(ethnicity_groups_imputed = if_else(ethnicity_groups_imputed |> is.na(), ethnicity_groups, ethnicity_groups_imputed)) |> 
      
      # Remove confounders of non interest
      tidybulk:::.resolve_complete_confounders_of_non_interest_df(dataset_id, assay_groups, disease_groups) |> 
      
      # Set intercept
      mutate(
        ethnicity_groups_imputed = fct_relevel(ethnicity_groups_imputed, "European"),
        assay_groups___altered = fct_relevel(assay_groups___altered, "10x Genomics 3"),
        disease_groups___altered = fct_relevel(disease_groups___altered, "Normal"),
        age_bin = fct_relevel(age_bin, "Senior_50"),
        age_decade = fct_relevel(age_decade, "5")
      ) 
    
  }
  
  flag_enriched_samples_within_immune_compartment = function(
    cell_tbl,
    sample_column = "sample_id",
    min_immune_cells = 50,
    min_reference_datasets = 3,
    entropy_z_threshold = -3,
    lineage_dominance = 0.9,
    lineage_reference_max = 0.5
  ) {
    
    # [HEURISTIC] Tissue groups treated as brain for the brain_nuclei rule
    brain_tissues = c(
      "cerebral lobes and cortical areas", "brainstem and cerebellar structures",
      "limbic and basal systems", "general brain and major structures"
    )
    
    # [HEURISTIC] Metadata text match for sorted populations. Lineage words must appear
    # within the first five words of the title, which is where sorted-population
    # datasets name themselves ("B cells from ..."), not in study descriptions
    # ("... tissue-resident T cells in ccRCC")
    sorted_title_regex = paste0(
      "^(\\S+\\s+){0,4}(B|T|NK|NKT|natural killer( T)?|plasma|macrophage|monocyte|",
      "myeloid|dendritic|mast|neutrophil) cells\\b|microglia"
    )
    
    # [DATA] suspension_type is optional; without it brain nuclei cannot be identified
    if (!"suspension_type" %in% colnames(cell_tbl)) {
      warning("suspension_type not found: the brain_nuclei rule is skipped")
      cell_tbl = cell_tbl |> mutate(suspension_type = NA_character_)
    }
    
    # [DATA] Immune cell counts per sample x cell type, computed in the database so only
    # a small table is collected. Non-immune cells are ignored entirely: the immune
    # fraction of a sample does not affect the immune composition being modelled.
    # Cell types are also grouped into lineages (T, B, monocytic, DC, NK, ...)
    # following L0 of caq_celltype_level_map.
    counts =
      cell_tbl |>
      filter(is_immune, cell_type_unified_ensemble != "non_immune") |>
      rename(sample = all_of(sample_column)) |>
      count(sample, dataset_id, tissue_groups, title, suspension_type, cell_type_unified_ensemble) |>
      collect() |>
      mutate(
        n = as.numeric(n),
        lineage = case_when(
          cell_type_unified_ensemble %in% c("b memory", "b naive", "plasma") ~ "b",
          cell_type_unified_ensemble %in% c("cd14 mono", "cd16 mono") ~ "monocytic",
          cell_type_unified_ensemble %in% c("cdc", "pdc") ~ "dc",
          grepl("^cd4|^cd8|^treg$|^tgd$|^mait$", cell_type_unified_ensemble) ~ "t",
          TRUE ~ cell_type_unified_ensemble
        )
      )
    
    # [MATH] Shannon entropy of the immune composition, H = -sum_k p_k log p_k (nats).
    # H = 0 when all immune cells are one type, H = log K when K types are equally
    # abundant; enrichment for a subset concentrates the composition and lowers H.
    # [STATISTIC] The plug-in estimate is biased downwards by about (K - 1) / (2N) for
    # N cells and K observed types; the Miller-Madow correction adds this term back,
    # so small datasets do not look artificially "enriched".
    shannon_entropy = function(n) {
      n = n[n > 0]
      p = n / sum(n)
      -sum(p * log(p)) + (length(n) - 1) / (2 * sum(n))
    }
    
    # [STATISTIC] Leave-one-dataset-out reference: each dataset is compared with the
    # other datasets of the same tissue only, so a strongly enriched dataset cannot
    # shift its own reference. Tissues differ hugely in baseline immune composition
    # (brain ~100% microglia, blood diverse), hence no global reference.
    # [HEURISTIC] At least min_reference_datasets usable datasets are required,
    # otherwise the reference is NA and the data-driven rules do not fire.
    leave_one_out = function(x, usable, statistic) {
      map_dbl(seq_along(x), \(i) {
        reference = x[-i][usable[-i]]
        if (length(reference) < min_reference_datasets) NA_real_ else statistic(reference)
      })
    }
    
    # [MATH] Lineage shares per dataset x tissue, pooling cells across samples
    # (dataset x tissue is the unit because enrichment is a protocol choice).
    # [STATISTIC] Reference share of each lineage = leave-one-out median across the
    # tissue's other datasets (median: robust to other enriched datasets).
    # Keeps the dominant lineage of each dataset with its reference share.
    lineage_share =
      counts |>
      count(dataset_id, tissue_groups, lineage, wt = n) |>
      complete(nesting(dataset_id, tissue_groups), lineage, fill = list(n = 0)) |>
      group_by(dataset_id, tissue_groups) |>
      mutate(share = n / sum(n), usable = sum(n) >= min_immune_cells) |>
      group_by(tissue_groups, lineage) |>
      mutate(reference_share = leave_one_out(share, usable, median)) |>
      group_by(dataset_id, tissue_groups) |>
      slice_max(share, n = 1, with_ties = FALSE) |>
      ungroup() |>
      select(dataset_id, tissue_groups, top_lineage = lineage, top_lineage_share = share, reference_share)
    
    # [MATH] Pooled immune entropy per dataset x tissue.
    # [STATISTIC] Reference = leave-one-out median and MAD (scaled to be consistent
    # with the SD under normality) of the tissue's other datasets.
    # [HEURISTIC] min_immune_cells: datasets with fewer immune cells are neither
    # assessed nor used as reference. The MAD is floored at 0.15 nats so tissues
    # whose datasets all look alike (e.g. brain, MAD ~ 0) do not flag tiny deviations.
    datasets =
      counts |>
      group_by(dataset_id, tissue_groups) |>
      summarise(
        title = first(title),
        nuclei = any(suspension_type %in% "nucleus"),
        n_immune = sum(n),
        entropy = shannon_entropy(tapply(n, cell_type_unified_ensemble, sum)),
        .groups = "drop"
      ) |>
      group_by(tissue_groups) |>
      mutate(
        usable = n_immune >= min_immune_cells,
        reference_entropy = leave_one_out(entropy, usable, median),
        reference_mad = pmax(leave_one_out(entropy, usable, mad), 0.15)
      ) |>
      ungroup() |>
      left_join(lineage_share, by = c("dataset_id", "tissue_groups")) |>
      
      # Decision rules; a dataset is dropped if any rule fires
      # low_entropy   [STATISTIC] robust z-score z = (H - median_ref) / MAD_ref;
      #               [HEURISTIC] cutoff z <= -3, the conventional robust outlier bound
      # one_lineage   [HEURISTIC] one lineage >= 90% of immune cells where its typical
      #               share in the tissue is <= 50%. Catches lineage sorting even when
      #               the tissue reference is too variable for the z-score to fire, and
      #               spares tissue-intrinsic dominance (brain microglia)
      # sorted_title  [HEURISTIC] metadata evidence of sorting, independent of counts
      #               (catches partial sorts that leave the composition diverse)
      # brain_nuclei  [HEURISTIC] not enrichment: brain snRNA-seq immune compartment is
      #               ~99% microglia, so these samples carry no compositional
      #               information; dropping them avoids false macrophage zeros that
      #               removing microglia from brain only would create
      mutate(
        entropy_z = (entropy - reference_entropy) / reference_mad,
        low_entropy = usable & coalesce(entropy_z <= entropy_z_threshold, FALSE),
        one_lineage = usable & top_lineage_share >= lineage_dominance & coalesce(reference_share <= lineage_reference_max, FALSE),
        sorted_title = coalesce(grepl(sorted_title_regex, title, ignore.case = TRUE, perl = TRUE), FALSE),
        brain_nuclei = nuclei & tissue_groups %in% brain_tissues,
        drop = low_entropy | one_lineage | sorted_title | brain_nuclei
      )
    
    # [DATA] Every sample inherits the decision of its dataset x tissue. A sample
    # spanning several (e.g. sample_column = "donor_id") is flagged if any is flagged
    counts |>
      left_join(datasets |> select(dataset_id, tissue_groups, drop), by = c("dataset_id", "tissue_groups")) |>
      summarise(flag = any(drop), .by = sample) |>
      rename(!!sample_column := sample)
  }
  

  
  pipeline = 
    tt_initialise(
    store = result_directory,
    computing_resources = tt_controller_elastic_slurm(
      workers = 20,
      time_hours = 28,
      mem_gb_per_job = c(5, 10, 20, 50, 80, 100, 150, 200), 
      cpus_per_task = 28,
      crashes_max = 1,
      script_lines = "export LD_LIBRARY_PATH=$HOME/envs/openssl/lib:$LD_LIBRARY_PATH"
    ),
    # debug_step = "input_relative_sample",
    garbage_collection = 100,
    workspace_on_error = TRUE,
    error = "stop",
    packages = c(
      "tidytargets", "targets", "tarchetypes", "tidyverse", "glue", "qs2",
      "dplyr", "tidyr", "tibble", "purrr", "stringr", "readr", "forcats",
      "magrittr", "duckdb", "dbplyr", "tidybulk", "cellNexus",
      "cellxgene.census", "sccomp", "crew", "crew.cluster"
    )
  ) |>
    tt_global(
      check_rclone_installation, rclone_copy, rclone_read_csv,
      edit_covariates, create_input_cell_counts, flag_enriched_samples_within_immune_compartment,
      result_directory
    ) |>
    tt_single(
      ethnicity_imputed <- rclone_read_csv(
        "box:/Mangiola_ImmuneAtlas/reports/ning/data/All_pseudobulk_1_0_6_ethnicity_imputed_colData.csv"
      ) |>
        select(sample_id, ethnicity_groups, ethnicity_groups_imputed = finalEthnicity_groups) |>
        mutate(ethnicity_groups_imputed = ethnicity_groups_imputed |> str_replace("_imp$", "_imputed"))
    ) |>
    tt_single(
      caq_celltype_level_map <- rclone_read_csv(
        "box:/Mangiola_ImmuneAtlas/reannotation_consensus/caq_celltype_level_map.csv"
      )
    ) |>
    tt_single(
      input_relative_sample <- create_input_cell_counts(
        caq_celltype_level_map, ethnicity_imputed,
        count_by = "sample_id"
      ), 
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_100"))
    ) |>
    tt_single(
      input_relative_donor <- create_input_cell_counts(
        caq_celltype_level_map, ethnicity_imputed,
        count_by = "donor_id"
      ),
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_100"))
    ) |>
    tt_single(
      saved_input_relative <- {
        file_name <- glue("{result_directory}/cell_metadata_1_0_15_sccomp_input_counts.rds")
        saveRDS(input_relative_sample, file_name)
        rclone_copy(file_name, "box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
      },
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_100"))
    ) |>
    tt_split(
      formula_df <- tribble(
        ~ name, ~ formula_composition, ~ formula_variability, ~ counts,
        
        # continuous age
        "estimates_continuous_age",
        "~ 1 + age_days_scaled*sex + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_days_scaled*sex + ethnicity_groups_imputed | tissue_groups)",  
        "~ 1",
        input_relative_sample,
        
        # discrete
        "estimates_age_bins",
        "~ 1 + age_bin + sex + age_bin:sex + ethnicity_groups_imputed + assay_groups___altered +
          (1 | dataset_id___altered) +
          (1 + age_bin + sex + age_bin:sex + ethnicity_groups_imputed | tissue_groups)",
        "~ 1",
        input_relative_sample,
        
        # discrete decade
        "estimates_age_decade", 
        "~ 1 + age_decade + sex + age_decade:sex + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_decade + sex + age_decade:sex + ethnicity_groups_imputed | tissue_groups)",  
        "~ 1 ",
        input_relative_sample,
        
        # discrete decade
        "estimates_age_decade_plasma", 
        "~ 1 + age_decade + sex + age_decade:sex + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered) + 
          (1 + age_decade + sex + age_decade:sex + ethnicity_groups_imputed | tissue_groups)",  
        "~ 1",
        input_relative_sample ,
        
        # interaction decade, this model is to test that the signal really comes from each tissue and not learned
        # at the body level, because of lack of data
        "estimates_age_decade_interaction", 
        "~ 1 + age_decade*sex*tissue_groups + ethnicity_groups_imputed + assay_groups___altered + 
          (1 | dataset_id___altered)", 
        "~ 1",
        input_relative_sample,
        
        # discrete decade DONOR grouping
        "estimates_age_decade_donor_grouping",
        "~ 1 + age_decade + sex + age_decade:sex + ethnicity_groups_imputed + assay_groups___altered +
          (1 | dataset_id___altered) +
          (1 + age_decade + sex + age_decade:sex + ethnicity_groups_imputed | tissue_groups)",
        "~ 1",
        input_relative_donor,
        
        # interaction decade DONOR grouping
        "estimates_age_decade_interaction_donor_grouping",
        "~ 1 + age_decade*sex*tissue_groups + ethnicity_groups_imputed + assay_groups___altered +
          (1 | dataset_id___altered)",
        "~ 1", # + ",
        input_relative_donor,
        
        # discrete + interaction ethnicity sex
        "estimates_age_bins_sex_ethnicity_interaction",
        "~ 1 + age_bin + sex + age_bin:sex + ethnicity_groups_imputed * sex + assay_groups___altered +
          (1 | dataset_id___altered) +
          (1 + age_bin + sex + age_bin:sex + ethnicity_groups_imputed * sex | tissue_groups)",
        "~ 1 ",
        input_relative_sample,
        
        # continuous + discrete
        "estimates_continuous_age_plus_age_bins",
        "~ 1 + age_days_scaled + age_bin*sex + ethnicity_groups_imputed + assay_groups___altered +
          (1 | dataset_id___altered) +
          (1 + age_days_scaled + age_bin*sex + ethnicity_groups_imputed | tissue_groups)",
        "~ 1 ",
        input_relative_sample,
        
        # smooths
        # Sex-specific age curves per tissue, so the age effect is modulated by
        # sex exactly as `age_decade:sex` does in estimates_age_decade and the
        # two models answer the same question.
        #
        # `bs = "fs"` takes one factor only, so tissue and sex are folded into
        # `tissue_sex` (built in the `estimates` target below). This costs no
        # extra random-effect slots: the number of penalty blocks a factor
        # smooth needs depends only on `k`, not on how many levels the factor
        # has, so the six tissue x sex cells still take 3 blocks. Total 5 of 6
        # slots: 1 dataset + 1 tissue + 3 fs blocks.
        #
        # `sex` is deliberately absent from the tissue clause. The fs basis
        # spans the per-level constants of its factor exactly, so this smooth
        # already owns all six tissue x sex baselines; a `sex` term in the
        # tissue clause would compete for the same quantities and only the sum
        # would be identified. Tissue-specific sex baselines are still in the
        # model, they just live in the smooth.
        #
        # Ethnicity is `1 +` rather than `0 +` for the same reason: the
        # dummies of a `0 +` block sum to one per sample, so that block spans
        # the tissue intercept and each tissue's baseline moves by an
        # arbitrary amount decided by the two shrinkage SDs.
        "estimates_smooths",
        "~ 1 + sex +
          ethnicity_groups_imputed + 
          assay_groups___altered +
          s(age_days, k = 8) + 
          s(age_days, tissue_sex, bs = \"fs\", k = 8) +
          (0 + ethnicity_groups_imputed | tissue_groups) +
          (1 | dataset_id___altered)",
        
        "~ 1",
        input_relative_sample,
        
        # disease tissue specific
        "estimates_age_bins_disease",
        "~ 1 + age_bin + sex + age_bin:sex + ethnicity_groups_imputed + assay_groups___altered +
          (1 | dataset_id___altered) +
          (1 + age_bin + sex + age_bin:sex + ethnicity_groups_imputed | tissue_groups)",
        "~ 1",
        input_relative_sample
      ) |> 
        expand_grid(
          cell_type_level = glue("L{0:4}") |> as.character(), 
          drop_disease = TRUE, # c(TRUE, FALSE),
          immune_only = c(TRUE,FALSE)
        ) |> 
        
        # Keep NON immune for one model only
        filter(immune_only | name == "estimates_age_decade") |> 
        
        # Filter immune if needed
        mutate(counts = map2(counts, immune_only, ~ {
          if (.y) {
            .x |> filter(is_immune)
          } else {
            .x |>
              nest(data = -c(sample_id, tissue_groups)) |>
              filter(map_lgl(data, ~ any(.x$cell_type_unified_ensemble == "non_immune"))) |>
              mutate(
                non_immune_count = map_int(data, ~ sum(.x$n[.x$cell_type_unified_ensemble == "non_immune"])),
                total_count = map_int(data, ~ sum(.x$n)),
                proportion_non_immune = non_immune_count / total_count
              ) |>
              filter(proportion_non_immune > 2 / 3) |>
              filter(!tissue_groups %in% c(
                "bone marrow",
                "lymphatic system",
                "spleen",
                "thymus",
                "blood"
              )) |>
              unnest(data)
          }
        })) |> 
        
        # Drop Disease if needed
        mutate(
          counts = map2(counts, drop_disease, 
                        ~ {
                          if(.y) .x |> filter(disease_groups == "Normal")
                          else .x 
                        }),
          formula_composition = if_else(drop_disease, formula_composition |> str_remove_all("\\+ disease_groups___altered"), formula_composition),
          formula_variability = if_else(drop_disease, formula_variability |> str_remove_all("\\+ disease_groups___altered"), formula_variability)
        ) |> 
        mutate(local_file_name = glue("{name}___{cell_type_level}___disease_{!drop_disease}___immune_only_{immune_only}")) |> 
        
        # TEMPORARY
        # FILTER FOR JUST ONE MODEL
        filter(name %in% c("estimates_age_decade", "estimates_age_decade_interaction",  "estimates_age_decade_interaction_donor_grouping", "estimates_smooths")) |> 
        #filter(name %in% c("estimates_age_decade", "estimates_age_decade_plasma", "estimates_age_decade_interaction", "estimates_age_decade_donor_grouping", "estimates_age_decade_interaction_donor_grouping")) |> 
        filter(cell_type_level %in% c( "L3") | (cell_type_level == "L4" & name == "estimates_age_decade_plasma")) |> 
        group_split(local_file_name)
    )  |>
    tt_iterate(
      estimates <- {
        formula_df$counts[[1]] |>
          
          # One factor per tissue x sex cell, for the fs smooth in
          # estimates_smooths. `bs = "fs"` accepts a single factor, so this is
          # what makes its age curve sex specific. Unused by the other models.
          mutate(tissue_sex = interaction(tissue_groups, sex, drop = TRUE)) |>
          
          # With L0 I have to summarise further because I have counts already
          with_groups(
            c(sample_id,age_days, age_days_scaled,age_bin, age_decade, sex, disease_groups___altered,ethnicity_groups_imputed, assay_groups___altered, dataset_id___altered, tissue_groups, tissue_sex, all_of(formula_df$cell_type_level[[1]])),
            ~ .x |> summarise(n = sum(n))
          ) |>
          #
          sccomp_estimate(
            formula_composition = formula_df$formula_composition[[1]] |> as.formula(),
            formula_variability = formula_df$formula_variability[[1]] |> as.formula(),        # Differential variability
            sample = "sample_id",
            cell_group = formula_df$cell_type_level[[1]], # A level of the hierarchy
            abundance = "n",
            cores = parallelly::availableCores(),
            mcmc_seed = 42,
            verbose = T,
            bimodal_mean_variability_association = TRUE,
            # prior_mean = list(intercept = c(0, 0.8), coefficients = c(0, 3)),
            # prior_overdispersion_mean_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40)),
            output_directory = "/vast/scratch/users/shen.m/sccomp_draws_files",
            
            max_sampling_iterations = 1000,
            warmup_samples = 400,
            
            inference_method =   "hmc",
            refresh = 1
          )
      },
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_150"))


    ) |>
    tt_iterate(
      saved_and_tranferred <- {
        fit <- estimates
        name <- formula_df$local_file_name[[1]]
        fit_path <- file.path(result_directory, paste0(name, "_FIT_FOR_PORTABILITY.rds"))
        out_path <- file.path(result_directory, paste0(name, ".rds"))
        
        attr(fit, "fit")$save_object(file = fit_path)
        rclone_copy(
          fit_path,
          "box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_estimates_1_0_15/"
        )
        attr(fit, "fit") <- readRDS(fit_path)
        
        fit |> sccomp_test() |> saveRDS(out_path)
        rclone_copy(
          out_path,
          "box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_on_cellNexus_1_0_15_add_smooths/"
        )
      },
      resources = tar_resources(crew = tar_resources_crew(controller = "elastic_20"))
    ) 
  
  pipeline |> tt_evaluate()

  
})



system("/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone copy /vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_bins___L3___disease_TRUE___immune_only_TRUE.rds box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/sccomp_on_cellNexus_1_0_15_hca_2025/")

pipeline |> tidytargets::tt_explore("formula_df")



library(targets)

input_relative_sample = tar_read(input_relative_sample, store = "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets")
input_relative_donor = tar_read(input_relative_donor, store = "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets")


tar_workspace(input_relative_sample, 
              script = "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets.R", 
              store = "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets"
)

tar_workspace("estimates_4bdedc368a84e9da", store = "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets", script = "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets.R")
tar_workspace("input_relative_sample", store = "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets", script = "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets.R")


library(tidyverse)
library(sccomp)
library(magrittr)
library(glue)
library(forcats)
library(stringr)
library(arrow)
library(dplyr)
library(duckdb)

tar_meta(starts_with("estimates_"), store = "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets")

# Age proportion prediciton
estimate_age_decade = readRDS("sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_decade_interaction_donor_grouping___L3___disease_FALSE___immune_only_TRUE.rds")

estimate_age_decade |> 
  sccomp_predict(
    formula_composition = ~ 1 + age_decade*sex + (1 + age_decade*sex | tissue_groups), 
    number_of_draws = 100,
    summary_instead_of_draws = TRUE
  ) |> 
  # |> 
  #   mutate(age_bin = factor(
  #     age_bin, 
  #     c("Infancy", "Childhood", "Adolescence", "Young Adulthood", "Middle Age", "Senior"),
  #     ordered = TRUE
  #   )) |> 
  #   mutate(age_bin_numeric = age_bin |> as.integer())  |> 
  saveRDS("/vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/prediction_age_decade.rds")

prediction_proportion = readRDS("sccomp_on_cellNexus_1_0_15_hca_2025/")


system("/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone copy /vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/prediction_age_bins.rds box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")

tar_read(formula_df, store = "/vast/scratch/users/shen.m/sccomp_on_cellNexus_hca_2025/_targets")


# For Hong
estimate_age_decade = readRDS("sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_decade_interaction___L3___disease_FALSE___immune_only_TRUE.rds")
estimate_age_bins = estimate_age_bins |> dplyr::select(-count_data)
attr(estimate_age_bins, "fit") = NULL
estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_bins_effect_tibble_only.rds")
system("/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone copy /vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_bins_effect_tibble_only.rds box:/immune_map_disease/")

# Save fit
library(magrittr)
estimate_age_bins = readRDS("/vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_bins.rds")
estimate_age_bins |> attr("fit") %$% save_object(file = "/vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_bins_FIT_FOR_PORTABILITY.rds") 
system("/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone copy /vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_bins_FIT_FOR_PORTABILITY.rds box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")
estimate_age_bins |> attr("fit") = readRDS("/vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_bins_FIT_FOR_PORTABILITY.rds")
estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_bins.rds")
system("/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone copy /vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_bins.rds box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")


# estimate_age_bins |> saveRDS("/vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/estimates_age_bins.rds")
system("/stornext/System/data/software/rhel/9/base/datamover/rclone/1.67.0/bin/rclone copy /vast/projects/mangiola_immune_map/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_15_hca_2025/21_11_2024_sccomp_archive_before_factor_ordering/estimates_age_bins.rds box:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")

# Benchmark
tic()
estimate_age_bins |> sccomp_test(contrasts = c(  "respiratory system" = "sexmale + `sexmale___respiratory system`",
                                                 "blood" = "sexmale + `sexmale___blood`"))
toc()


estimate_age_bins |> 
  sccomp_test()



# annotate plasma
out_dir <- file.path(getwd(), "./", "plasma_h5se")


build_plasma_sce <- function(cache_dir,
                             subset_n = 100,
                             min_features = 8000,
                             seed = NULL) {
  metadata <- cellNexus::get_metadata(
    cloud_metadata = get_metadata_url(c("hca_2024", "hca_2025"))
  ) 
  #|> 
  # left_join(get_census_table())
  
  plasma_metadata <- metadata |>
    dplyr::filter(cell_type_unified_ensemble == "plasma") |>
    dplyr::filter(
      empty_droplet == FALSE,
      alive == TRUE,
      scDblFinder.class != "doublet",
      feature_count >= min_features
    ) |>
    dplyr::collect()
  
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(subset_n) && nrow(plasma_metadata) > subset_n) {
    plasma_metadata <- plasma_metadata |>
      head(subset_n)
  }
  
  sce <- cellNexus::get_single_cell_experiment(
    plasma_metadata,
    cache_directory = cache_dir
  )
  
  sce
}


compute_singscore_delayed <- function(
    sce, assay_name = NULL, upSet = singscore::tgfb_gs_up, 
    downSet = singscore::tgfb_gs_dn, workers = 1, output_file = NULL
) {
  assay_names <- SummarizedExperiment::assayNames(sce)
  if (is.null(assay_name)) {
    assay_name <- if ("logcounts" %in% assay_names) "logcounts" else assay_names[[1]]
  }
  message("Ranking genes using assay: ", assay_name)
  expr_mat <- SummarizedExperiment::assay(sce, assay_name)
  feature_ids <- rownames(expr_mat)
  
  # Debug information
  message("Expression matrix dimensions: ", nrow(expr_mat), " x ", ncol(expr_mat))
  message("Expression matrix rownames length: ", length(rownames(expr_mat)))
  message("Expression matrix class: ", class(expr_mat))
  # Choose ID space (prefer best overlap with matrix rownames)
  up_sym <- upSet; down_sym <- downSet
  # Map symbols -> ENSG
  up_ens <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = up_sym, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first")
  down_ens <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = down_sym, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first")
  
  # Use delayed ranking with optional output file
  ranked <- singscore::rankGenes(
    expr_mat, 
    workers = workers 
    
    # TEMPORARY
    # , stableGenes = up_ens    
  )
  
  # Convert character vectors to GeneSet objects
  library(GSEABase)
  up_geneset <- GeneSet(up_ens, setName = "upSet")
  down_geneset <- GeneSet(down_ens, setName = "downSet")
  
  # Use regular simpleScore for now since the delayed version has issues
  # Convert DelayedMatrix to matrix for simpleScore
  # Set the stable attribute properly
  scoredf_delayed <- singscore::simpleScore(
    ranked, 
    upSet = up_geneset, 
    downSet = down_geneset
  )
  
  # Add original cell
  scoredf_delayed$cell_id = sce |> colData() |> _[,"original_cell_"]
  scoredf_delayed$dataset_id = sce |> colData() |> _[,"dataset_id"]
  
  # Return the delayed results and the ranked data
  return(scoredf_delayed)
}


# Build from cache (default)
sce <- build_plasma_sce(cache_dir = "/vast/projects/cellxgene_curated/cellNexus", subset_n = Inf, min_features = 8000) 

library(tidySummarizedExperiment)
# BiocManager::install("stemangiola/singscore@delayedArray")
library(singscore)
# Rename features using

# from sce to singscore execution
long_lived_genes <- c("SDC1", "CD38", "TNFRSF17", "TNFRSF13B", "CD28", "CXCR4", "CD69", "ITGA4", "ITGB1", "ITGB7", "ITGAE", "PRDM1", "XBP1", "MZB1", "DERL3", "HSPA5", "BCL2", "MCL1")
short_lived_genes <- c("CD19", "PTPRC", "HLA-DRA", "HLA-DRB1", "CD27", "CXCR3", "S1PR1", "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "PAX5", "BCL6", "BACH2", "AICDA", "CCR9", "CCR10")

long_lived_genes <- c("SDC1", "TNFRSF17", "CXCR4", "PRDM1", "XBP1", "MZB1", "MCL1")
short_lived_genes <- c("CD19", "PTPRC", "HLA-DRA", "MKI67", "PAX5", "BCL6", "BACH2")

sce_scored_delayed <- compute_singscore_delayed(
  sce, 
  assay_name = "counts", 
  upSet =long_lived_genes  , 
  downSet = short_lived_genes, 
  workers = 20,
  output_file = NULL  # Disable output file for now to avoid HDF5 conflicts
)

sce_scored_delayed |> saveRDS("sce_scored_delayed_plasma.rds")

sce_scored_delayed |>
  as_tibble(rownames="cell") |> 
  left_join(colData(sce) |> as_tibble(rownames="cell")) |>
  ggplot(aes(x=TotalScore)) +
  geom_density(aes(color = dataset_id)) +
  guides(color = "none") +
  ylim(c(0, 30)) + 
  theme_minimal() 

plot(sce_scored_delayed$TotalScore, sce_scored_delayed_down$TotalScore)

# sce_scored |>
#colData() |>
#saveRDS(file = "dev/plasma_h5se_scored.rds")

