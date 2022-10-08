# Import libraries
import streamlit as st
from streamlit_option_menu import option_menu
import pandas as pd
import leafmap.foliumap as leafmap
from stmol import showmol
import py3Dmol
from Bio import SeqIO, Entrez
from Bio.SeqUtils import GC, molecular_weight as mw,MeltingTemp as mt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
import matplotlib.pyplot as plt
from http.client import IncompleteRead
# Page config
st.set_page_config(layout='wide')
# Menu bar
with st.sidebar:
    select=option_menu(menu_title='Menu',options=['Overview','Genomics & Proteomics','References'],icons=['book-half','infinity','journal-text','file-earmark-person'],default_index=0,styles={'container':{'background-color':'#FFFFFF'},'nav-link':{'font-size':'14px','text-align':'left'},'nav-link-selected':{'background-color':'#0557FD'}})
#############################################################################################################################################The Overview
############################################################################################################################################
if select=='Overview':
    st.sidebar.markdown('***')
    #About
    st.sidebar.subheader('About')
    st.sidebar.info('''
    This app is made using :\n 
    Streamlit - Python - Biopython - Pandas- Py3dmol - Matplotlib - Leafmap\n\n
    It is maintained by [Yassine OUCHEN](https://github.com/ouchen-bio/Anammox-Bacteria-Database)
    Email: yassine.ouchen2021@gmail.com
    LinkedIn: [Click Here](https://www.linkedin.com/in/yassine-ouchen-be)
    ''')
    # Title
    cl1,cl2,col3=st.columns([1,2.4,1])
    with cl2:
        st.title('The World Of Anammox')   
    st.markdown('***')
    #Introduction
    st.header('Introduction')
    st.markdown('')
    st.write('Over the last decades, the excess release of nitrogen into water has increased drastically, and nitrogen should be immediately removed from effluents before their discharge into the environment. Conventionally, biological nitrogen elimination from wastewater has been achieved through full **nitrification-denitrification** process. During nitrification phase, ammonium is oxidized to nitrite by ammonium oxidizing bacteria and then to nitrate by nitrite oxidizing bacteria. During denitrification, nitrates are denitrified to dinitrogen gas under anoxic conditions in the presence of organic carbon compounds. Recently, a new process named **Anammox** (anaerobic ammonium oxidizing) has become an innovative technological advancement for industrial wastewater treatment plants. This process was first discovered in a denitrifying bioreactor in the early nineties (*Mulder et al. 1995*) and has great potential for the elimination of ammonia nitrogen from wastewater. The responsible bacteria directly convert ammonium and nitric oxide into dinitrogen gas under anoxic conditions. The application of anammox process is considered to be an efficient, eco-friendly and low-cost alternative to conventional biological nitrogen removal processes (*Jetten et al. 1997*).')
    # Anammox Bacteria
    st.markdown('')
    st.header('Anammox Bacteria')
    st.markdown('')
    st.subheader('1- Anammox Species Diversity')
    st.write("The bacteria that perform the anammox process are genera that belong to the bacterial phylum **Planctomycetes**. Until now, eighteen species have been confirmed. Because none of these have been successfully isolated or grown in a classical pure culture, all have been described as  **‘Candidatus’** (Ca) according to recommendations of the International Committee on Bacterial Nomenclature (*Sneath 1990*; *Murray and Stackebrandt 1995*). The eighteen anammox species are divided over sex genera as follow:")
    st.write("1. **Kuenenia**, one species: Ca. Kuenenia stuttgartiensis (*Strous et al. 2006*).")
    st.write("2. **Brocadia**, five species: Ca. Brocadia anammoxidans (*Strous et al. 1999a*), Ca. Brocadia sinica (*Oshiki et al. 2011*), Ca. Brocadia fulgida (*Kartal et al. 2008*), Ca. Brocadia caroliniensis (*Magrí et al. 2012*), Ca. Brocadia sapporoensis (*Narita et al. 2017*).")
    st.write("3. **Jettenia**, three species: Ca. Jettenia asiatica (*Quan et al. 2008*), Ca. Jettenia caeni (*Ali et al. 2015*), Ca. Jettenia ecosi (*Botchkova et al. 2018*)")
    st.write("4. **Anammoxoglobus**, one species: Ca. Anammoxoglobus propionicus (*Kartal et al. 2007a*).")
    st.write("5. **Scalindua**, seven species: Ca. Scalindua brodae, Ca. Scalindua wagneri (*Schmid et al. 2003*), Ca. Scalindua sorokinii (*Woebken et al. 2008*), Ca. Scalindua marina (*Brandsma et al. 2011*), Ca. Scalindua profunda (*van de Vossenberg et al. 2013*), Ca. Scalindua arabica (*Woebken et al. 2008*), Ca. Scalindua Japonica (*Oshiki et al. 2017*).")
    st.write("6. **Anammoxomicrobium**, one species: Ca. Anammoxomicrobium moscowii (*Khramenkov et al. 2013*).")
    st.markdown('')
    st.subheader('2- Anammox Cell Biology')
    st.markdown('')
    st.write('**2-1 Cell Structure**')
    st.markdown('')
    st.write('Anammox bacterium cell is a highly compact sphere (cocci) with a diametre of about 0.8 mm when observed under a transmission electron microscope (TEM) (*Strous et al. 1999a* ; *Fuerst 2005*; *Fuerst & Sagulenko 2011*). The cell consists of a cell wall devoid of peptidogycan, a cytoplasmic membrane and an intracytoplasmic membrane. The cytoplasm is devided into three essential compartments: riboplasm containing the nucleoid, paryphoplasm separated from the riboplasm by an intracytoplasmic membrane and anammoxosome containing tubule-like structures (Fig. 1) (*Lindsay et al. 2001*; *Van Niftrik et al. 2008b*; *Van Niftrik & Jetten 2012*). This last compartment is specific to anammox bacteria.')
    cl1,cl2,cl3=st.columns([1,8,1])
    with cl2:
        st.markdown('***')
        st.image("data/figures/celanamx.jpg")
        st.markdown("***")
    cl1,clt,cl3=st.columns([1,1.5,1])
    with clt:
        st.write('**Figure 1** Cell structure of anammox bacteria.')
    st.markdown("")
    st.write('**2-2 Anammoxosome: The power plant of anammox bacteria**')
    st.markdown('')
    st.write('Anammoxosome is the largest compartment of anammox bacterium and harbours the responsible enzymes of its metabolism. This organelle is separated from the rest of the cell by a single bilayer membrane. The anammoxosome membrane possesses a unique type of lipids, named **ladderane**, with sequential structure of either five linearly condensed cyclobutanes or three cyclobutanes and one cyclohexane arranged like a "staircase" at the end of the hydrocarbon chains (*Kuypers et al. 2003, 2005*; *Damsté et al. 2002, 2005*; *Hopmans et al. 2006*). So far, four main classes of ladderane lipids have been identified in all anammox bacteria species: ladderane with an ester group linked to a heptyl, ladderane with an ether group bound to a glycerol molecule, ladderane with a combination of both ester-linked and ether-linked systems and a mixed ladderane/n-alkyl glycerol diether (*Shinnighe Damstré et al. 2004; Rattray et al. 2008*). Due to the structural properties of these lipids, the anammoxosome membrane can be function as an energy generator by maintaining a proton motive force and as a barrier against toxic intermediates of the anommox process (*Van Niftrik et al. 2004; Rattray et al. 2008*).')
    st.markdown('')
    #Nitrogen biotransformation
    st.header("Anammox reaction and its key enzymes")
    st.markdown("")
    st.write("Anaerobic ammonium-oxidizing bacteria are a group of chemolithoautotroph that oxidize ammonium in the presence of nitrite as an electron acceptor and produce dinitrogen gas as the end product under anaerobic conditions (Eq.1). It's well known that biological oxidation of ammonium usually occurs in aerobic or low-oxygen systems, so how these microorganisms are capable of performing such reaction anaerobically?. Undoubtedly, the process would proceed via intermediates and partial reactions. ")
    cl1,clt,cl3=st.columns([1,1,1])
    with clt:
        st.image("data/figures/main_equation.jpg",width=300)
    st.markdown('')
    st.write("At the beginning of anammox discovery, *Van de Graaf et al.(1997)* proposed a three-step process with hydroxylamine and hydrazine as key intermediates. The first step involved reduction of nitrite to hydroxylamine by hydroxylamine oxidoreductase, an enzyme that is present in aerobic ammonium-oxidizing bacteria and catalyzes the oxidation of hydroxylamine to nitrite (*Hooper et al. 1997*). The next step comprised the formation of hydrazine by condensation of hydroxylamine and ammonium. The last step was oxidation of hydrazine to generate dinitrogen gas. Later, genomics and proteomics analysis of K. stuttgartiensis's whole-cell revealed the presence of genes coding for nitric oxide oxidoreductase, which suggested nitric oxide (NO) to be an intermediate instead of hydroxylamine (NH2OH) (*Strous et al. 2006*). Based on these analysis, dinitrogen formation would be the result of the four-electron oxidation of hydrazine, catalyzed by an HAO-like enzyme, called hydrazine dehydrogenase (HDH) (Eq. 4). These four electrons would then lead to nitrite reduction by No-generating enzyme(s) (Eq. 2) and hydrazine production by hydrazine synthase (HZS) (Eq. 3).")
    cl1,clt,col3=st.columns([1,2.2,1])
    with clt:
        st.image('data/figures/reactions.jpg',width=500)
    st.markdown('')
    st.write('**>>> Hydrazine synthase (HZS)**') 
    st.write('HZS is a membre of the hydroxylamine oxidoreductase-like enzymes family that catalyzes the formation of hydrazine by condensation of ammonium and NO. The crystal structure of this multicomplex enzyme shows a dimer of heterotrimers αβγ where each of which has c-type cytochromes containing multiple active sites. According to the crystal structure, hydrazine syntesis would involve a reduction of nitric oxide to hydroxylamine at the active site of the γ-subnit, followed by a condensation of ammonium with hydroxylamine to produce hydrazine at the active site of the α-subnit (*Andreas Dielt et al. 2015*).')
    st.markdown('')
    st.write("**>>> No-generating enzymes**")
    st.write('In order to respond to different environmental conditions, anammox bacteria have developed different types of enzymes for NO production (*Kartal et al. 2013*; *Kartal & Keltjens 2016*). In K. stuttgartiensis a nitrite-reducing octaheme c-type cytochrome gene is highly expressed and was proposed to be involved in NO synthesis (*Kartal et al. 2013*; *Dietl et al., 2019*; *Ferousi et al. 2021*). The same organism contains a hardly transcribed gene that codes for another type of nitrite reductase called **NirS** (cytochrome cd1 nitrite reductase)(*Kartal et al. 2013*). In Scalindua profunda, on the other hand, both NirS and an octaheme c-type cytochrome are highly expressed, leading to conclude that this microorganism could employ two or more types of NO-producing enzymes (*van de Vossenberg et al. 2013*). In Jettenia caeni a copper-based **NirK-type nitrite reductase** was proposed to be involved in nitrite reduction for NO generation (*Hira et al. 2012*). Moreover, a recent study has expanded the list of No-producing enzymes in different anammox genra. This study revealed that the tubule-like structures of anammoxosome harbours a new different enzyme, named **nitrite oxidoreductase (NXR)** (an enzyme that is found in aerobic nitrite oxidizers ), which plays a key role in reduction of nitrite to nitric oxide (*Chicano et al. 2021*).')
    st.markdown('')
    st.write("**>>> Hydrazine dehydrogenase (HDH)**")
    st.write('HDH is a multiheme c-type cytochrome protein that oxidizes hydrazine, forming dinitrogen. HDH was first identified from the anammox species Kuenenia stuttgartiensis as the gene product of kustc0694 (*Kartal et al. 2011b*). This HAO-like enzyme occurs as an homotrimeric protein in which each monomer possesses eight c-type hemes. The heme 4 has been identified as the catalytic site for hydrazine oxidation which involves two main steps. The first step is a removal of 2 electrons from hydrazine (N2H4) to make Diimide (N2H2). The second step involves a removal of another 2 electons from Diimide to generate nitrogen gas (N2) (*Akram et al. 2019*).') 
    st.write('The following web-integrated application enables 3D visualization of anammox key enzymes. This app is powred by Py3Dmol.')
    #proteins 3d structure
    st.markdown("***")
    cly,clx=st.columns([2.62,.75])
    with clx:
        prot_name={'Hydrazine synthase':'5C2V','Hydrazine dehydrogenase':'6HIF','Cd1 nitrite reductase Nirs':'6TSI','Nitrite oxidoreductase':'7b04'}
        protein=st.selectbox('Select Enzyme :',prot_name)
        st.markdown('***')
        style=st.radio('Select Style :',['cartoon','line','sphere','stick'])
        st.markdown('***')
        b_color=st.color_picker('Set Background Color :','#F0F2F6')
        st.markdown('***')
    with cly:
        xyzview1=py3Dmol.view(query='pdb:'+prot_name.get(protein))
        xyzview1.setStyle({style:{'color':'spectrum'}})
        xyzview1.setBackgroundColor(b_color)
        showmol(xyzview1,height=500,width=800)
    st.markdown('***')
    #Physiology
    st.markdown('')
    st.header('Physiology')
    st.markdown('')
    st.write("**1- Physiological Characteristics**")
    st.markdown('')
    st.write('**Data Frame 1** Physiological characteristics of anammox bacteria (**n.d**: No data).')
    df_Ch=pd.read_csv('data/csv_files/anmx_phys.csv')
    df_Ch.set_index('physiological characteristics',inplace=True)
    filter_list=st.radio('',['all species','by species'],horizontal=True)
    if filter_list=='all species':
        st.dataframe(df_Ch)
    else:
        list_of_species=st.multiselect('',['Ca. Kuenenia Stuttgariensis','Ca. Brocadia sinica','Ca. Brocadia caroliniensis','Ca. Brocadia fulgida','Ca. Brocadia anammoxidans','Ca. Brocadia sapporoensis','Ca. Jettenia caeni','Ca. Jettenia ecosi','Ca. Scalindua sp','Ca. Scalindua profunda'])
        filtred_phy=df_Ch[list_of_species]
        st.dataframe(filtred_phy)
    st.markdown('')
    #Nitrogen removal chart
    st.markdown('')
    st.write("**2- Nitrogen Removal Rate**")
    df_k=pd.read_csv('data/csv_files/anmx_kenitics.csv')
    x=df_k['Time']
    y1=df_k['N-NH4 in']
    y2=df_k['N-NO2 in']
    y3=df_k['N-NH4 out']
    y4=df_k['N-NO3 out']
    y5=df_k['N-NO2 out']
    z1=df_k['N Load']
    z2=df_k['N Removal']
    z3=df_k['N Removal per']
    st.set_option('deprecation.showPyplotGlobalUse', False)
    c1,c2,c3=st.columns([6,.2,6.6])
    plt.rcParams['font.size']='11'
    def NH_removal():
        plt.plot(x,y1,marker='.',label='N-NH4 in')
        plt.plot(x,y2,marker='.',label='N-N02 in',color='r')
        plt.plot(x,y3,marker='.',label='N-NH4 out',color='b')
        plt.plot(x,y4,marker='.',label='N-N03 out',color='black')
        plt.plot(x,y5,marker='.',label='N-N02 out',color='purple')
        plt.xlabel('Time [days]')
        plt.ylabel('N concentration [g m−3 day−1]')
        plt.legend()
        plt.grid()
        plt.title('(a)')
        return st.pyplot()
    def N_removal():
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(x,z1,marker='.',label='N Load',color='brown')
        ax.plot(x,z2,marker='.',label='N Removal Rate',color='green')
        ax2=ax.twinx()
        ax2.plot(x,z3,marker='.',label='N Removal %',color='b')
        ax.set_xlabel('Time [days]')
        ax.set_ylabel('N Load and Removal Rate\n[kg m−3 day−1]')
        ax2.set_ylabel('N Removal Rate [%]')
        plt.title('(b)')
        ax.set_ylim(0,0.40)
        ax.legend(loc=0)
        ax2.legend(loc=0)
        ax.grid()
        return st.pyplot()
    with c1:
        NH_removal()
    with c3:
        N_removal()
    st.write('**Figure 2** An example of nitrogen removal performance of anammox bacteria in a SBR reactor during 85 days; **(a)** ammonia, nitrite and nitrate concentrations, **(b)** total nitrogen loads and their removal rates.')
    #Enrichement and identification of anammox bacteria
    st.header("Enrichement and identification of anammox bacteria")
    st.markdown("")
    st.subheader("1- Samples Collection")
    st.markdown("")
    st.write('Albeit the first anammox species (*B. Anammoxidans*) was discovred in a man-made system (wastewater), anammox bacteria can be found in any environment including soil, groundwater, freshwater, marine sediments, lakes, estuaries, oxygen minimum zones and continental shelves in the oceans, polar regions, hot springs, and deep-sea hydrothermal vents (*Op den Camp et al. 2006; Penton et al. 2006; Schmid et al. 2007; Jetten et al. 2009; Humbert et al. 2010*).')
    dfs=pd.read_csv('data/csv_files/anmx_source.csv')
    dfs.set_index('Species',inplace=True)
    st.markdown("")
    st.write("**Data Frame 2** Anammox bacteria species and their environmental sources.")
    st.dataframe(dfs)
    #interactive map
    in_csv='data/csv_files/amx_distribution.csv'
    Map=leafmap.Map()
    Map.add_points_from_xy(in_csv,x='longitude',y='latitude',color_column='system',icon_names=['leaf','recycle'],add_legend=True,spin=True)
    #Map.add_xy_data(in_csv,x='Longitude',y='Latitude')
    st.markdown('***')
    Map.to_streamlit()
    st.write('**Figure 3**: Interactive map presenting the global distribution of anammox sampling sites.')
    st.markdown('')
    #Enrichement of anammox bacteria
    st.subheader("2- Enrichement Of Anammox Bacteria")
    st.markdown('')
    st.write('**2-1- Inocula selection**')
    st.markdown("")
    st.write('Due to the extremely slow metabolism and long doubling time of anammox bacteria, different types of inoculum have been used for rapid start-up period  of the anammox process. These include denitrifying sludge (*van de Graaf et al. 1996*), anammox biomass (*Trigo et al. 2006; Casagrande et al. 2013*), nitrifying sludge (*van der Star et al. 2007*), marine sediments (*van de Vossenberg et al. 2008; Kindaichi et al. 2011*), activated sludge ( *Jin et al. 2011; Wang et al. 2011; Li et al. 2012* ) and mature granular anammox sludge (*Ni et al. 2011*).')
    st.markdown("")
    st.write('**2-2- Bioreactors**')
    st.markdown("")
    st.write('*>>> Sequencing batch reactor (SBR)*')
    st.markdown("")
    st.write('SBR is a type of activated sludge process for wastewater treatment that involves a sequence of steps that are performed in a single reactor tank. SBR has been reported to be an efficient system for anammox enrichment due to its unique ability to distribute the biomass homogeneously in the reactor. Because it has been stated that anammox bacteria growth is adversely affected by high concentration of nitrite, the enrichment in SBR provides a homogenous mixture of biomass and substrate in the reactor tank leading to the prevention of nitrite inhibition and therefore the increase of anammox activity (*Strous et al. 1998; Dapena-Mora et al. 2004; Arrojo & Mosquera-Corral 2006*).')
    st.markdown("")
    cl1,cl2,cl3=st.columns([1,1.2,1])
    with cl2:
        st.image("data/figures/Anammox_sbr.jpg",width=310)
        st.markdown("***")
    cl1,clt,cl3=st.columns([1,.8,1])
    with clt:
        st.write("**Figure 4** Anammox SBR reactor at Radboud university, Netherlands.")
    st.markdown("")
    st.write('*>>> Membrane bioreactor (MBR)*')
    st.markdown("")
    st.write(" MBR is another type of wastewater treatment reactor in which anammox bacteria are contained within the bioreactor by the application of a membrane system. Based on the location of the membrane MBRs are classified into submerged MBRs and external MBRs. Unlike SBR, MBR has less sludge production and a better biomass retention due to its impermeable membrane that acts as a microfilter. However, the higher cost for membrane material and the loss of membrane permeability limit the application of this type of bioreactors (*Trigo et al. 2006; Jiang et al. 2013*).")
    st.markdown("")
    st.write('*>>> Up-flow anaerobic biofilter (UAF)*')
    st.markdown('')
    st.write('UAF is an effective technology that has been widely employed for the treatment of a variety of wastewaters. The process involves the flow of wastewater through or over biomass growing on a packing media inside the bioreactor. The packing media provides a uniform flow through the bioreactor which improves the contact between substrates and the biomass, and promotes biomass retention of the anammox bacteria, resulting in a higher nitrogen removal rate (*Jin et al. 2008; Chen et al. 2010*).')
    st.markdown('')
    st.subheader('3- Anammox Bacteria Identification')
    st.markdown('')
    st.write('**3-1 Ladderane lipids measurement**')
    st.markdown('')
    st.write('In contrast to all other known microorganisms, anammox bacteria contain ladderane lipids in the cell membrane surrounding the anammoxosome. Because these special lipids are only found in anammox bacteria, they can be used as biomarkers to infer the presence of anammox bacteria in different environmental samples. The following steps present the protocol of ladderane lipids analysis as described by *Rattray et al. (2008)*:')
    st.markdown('')
    st.write('1. Samples of enrichment are collected from the bioreactor and anammox cells are obtained by precipitation using centrifugation.')
    st.write('2. The anammox cells are freeze-dried after removing and discarding the supernatant.')
    st.write('3. The dried biomass is extracted five times using a mixture of dichloromethane / methanol, and then the total lipid is obtained after removing the mixture by a rotatory evaporator.')
    st.write('4. The extracts are methylated with boron tri-fluoride / methanol, and then fractionated by a small silica column using ethyl acetate as the effluent. ')
    st.write('5. The final fraction is silylated with BSTFA (N,O-bis-(trimetylsilyl)trifluoroacetamide) in pyridine at 60 °C for 15 min.')
    st.write('6. A known amount of 6, 6-d2-3-methyle-icosane is added after silylation, and the derivatized extract is analyzed using gas chromatography - mass spectrometry (GC/MS).')
    st.write('**3-2 Fluorescence In Situ Hybridization**')
    st.markdown('')
    st.write('Fluorescence in situ hybridization, also known as FISH is a biomolecular technique that is used for the detection and quantification of nucleic acids in their cellular environment. The technique relies on the hybridization of a small labelled oligonucleotides sequences called probes to their complementary DNA target. FISH has been widely applied for quantifying and identifying anammox bacteria in different environmental samples (*Schmid et al. 2000; 2005*). Probes that target the 16S rRNA genes are the most commonly used in FISH detection of anammox bacteria, however probes for 23S rRNA gene and ISR (Intergenic Spacer Region) have also been applied (*Schmid et al. 2000 ; 2001*).')
    st.markdown('')
    st.write('**3-3 PCR-based method**')
    st.markdown('')
    st.write('Polymerase chain reaction is a molecular method that permits the amplification of a specific region of a given nucleic acid in order to detect and study it. PCR amplification is a three-step process that is carried out in repeated cycles: Denaturation of DNA double strands, annealing of the primers to complementary sequences and primers extension to synthesis new DNA strands. This biomolecular technique has been widely used to infer the presence of anammox bacteria in a wide range of samples. Most of anammox species have been detected and identified by using specific primers that target either the 16S rRNA genes or the functional genes (*Schmid et al. 2003; Schmid et al. 2008; Kuenen 2008; Lam et al. 2009; Jetten et al. 2009; Ni et al. 2010; Li et al. 2011b*).')
#############################################################################################################################################Genomics & Proteomics
############################################################################################################################################
elif select=='Genomics & Proteomics':
    # Read Data from NCBI
    def ANID_LTGN(AN):
        try:
            #read genebank File
            Entrez.email='anammox.bacteria18@gmail.com'
            handle=Entrez.efetch(db='nucleotide',id=AN, rettype='gb',retmode='text')
            seq_object=SeqIO.read(handle,"gb")
            #########Genomics#############
            st.markdown('***')
            st.header('Genomics')
            #DESCRIPTION
            definition=seq_object.description
            accession=seq_object.id
            Authors=seq_object.annotations['references'][0].authors
            Title=seq_object.annotations['references'][0].title
            Journal=seq_object.annotations['references'][0].journal
            complet_seq=seq_object.seq
            seq_lenght=len(complet_seq)
            gc_content=GC(complet_seq)
            #nucleotides number
            A_num=complet_seq.count("A")
            T_num=complet_seq.count("T")
            C_num=complet_seq.count("C")
            G_num=complet_seq.count("G")
            #nucleotides frequency
            A_freq=round((A_num/seq_lenght)*100,2)
            T_freq=round((T_num/seq_lenght)*100,2)
            C_freq=round((C_num/seq_lenght)*100,2)
            G_freq=round((G_num/seq_lenght)*100,2)
            st.subheader('**I- Genome (Total / Contig)**')
            st.write("**1- Description**")
            st.markdown('')
            st.write('**DEFINITION :**  ',definition)
            st.write('**ACCESSION :**   ',accession)
            st.write('**AUTHORS :**     ',Authors)
            st.write('**TITLE :**     ',Title)
            st.write('**JOURNAL :**     ',Journal)
            #features
            all_features=[feature.type for feature in seq_object.features]
            count_features=Counter(all_features)
            d={"Genome Length(bp)":seq_lenght,"GC Content(%)":[GC(complet_seq)],"gene(total)":count_features["gene"],"CDS(total)":count_features["CDS"],"tRNA":count_features['tRNA'],"rRNA":count_features["rRNA"],"ncRNA":count_features["ncRNA"],"tmRNA":count_features["tmRNA"]}
            info_table=pd.DataFrame(d)
            st.markdown('')
            ########Genome Info###########
            st.write("**2- Genome Info**")
            st.dataframe(info_table)
            st.markdown('')
            #Download full record
            full_seq=str(complet_seq)
            n=60
            splitted=[full_seq[i:i+n] for i in range(0,len(full_seq),n)]
            complet_seq_reslt='\n'.join(splitted)
            fasta_file=str('>'+accession+' '+definition+'\n'+complet_seq_reslt)
            st.write('Download Full Sequence: ')
            st.download_button(label='Download fasta',data=fasta_file,file_name='sequence.fasta',mime='text/plain')
            st.markdown('')
            ########Genome Info###########
            st.write("**3- Base Count**")
            nuc_num=dict([("Adenine",A_num/1000),("Thymine",T_num/1000),("Cytosine",C_num/1000),("Guanine",G_num/1000)])
            nuc_per=dict([("Adenine",A_freq),("Thymine",T_freq),("Cytosine",C_freq),("Guanine",G_freq)])
            nuc_filter=st.radio('',['Number','Percentage'],horizontal=True)
            st.set_option('deprecation.showPyplotGlobalUse', False)
            if nuc_filter=='Number':
                cx,cy=st.columns([3,1])
                with cx:
                    plt.rcParams['font.size']='8.5'
                    #plt.figure(figsize=(6,3.1))
                    dotn=plt.bar(nuc_num.keys(),nuc_num.values(),width=0.31,color=['green','blue','red','orange'])
                    for bar in dotn:
                        height=bar.get_height()
                        plt.annotate('{}'.format(height),
                             xy=(bar.get_x()+bar.get_width()/2,height),
                             va='bottom',
                             ha='center')
                    plt.xlabel('Nucleotides')
                    plt.ylabel('Number of nucleotides (Kbp)')
                    st.pyplot()
            else:
                cx,cy=st.columns([3,1])
                with cx:
                    plt.rcParams['font.size']='9'
                    dotp=plt.bar(nuc_per.keys(),nuc_per.values(),width=0.3,color=['green','blue','red','orange'])
                    for b in dotp:
                        height=b.get_height()
                        plt.annotate('{}'.format(height),
                             xy=(b.get_x()+b.get_width()/2,height),
                             va='bottom',
                             ha='center')
                    plt.xlabel('Nucleotides')
                    plt.ylabel('Percentage of nucleotides (%)')
                    st.pyplot()
            ############Genes##############
            st.write("**II- Genes**")
            st.markdown('')
            #Extract genes
            genes_features=[]
            for feature in seq_object.features:
                if feature.type=="gene":
                    genes_features.append(feature)
            #Empty Genes' parametres Lists
            genes_locustages=[]
            genes_names=[]
            genes_info=[]
            genes_sequences=[]
            genes_molecular_weight=[]
            st_melting_temp=[]
            all_g_fa=[]
            #loop through genes
            for genes in genes_features:
                if "locus_tag" in genes.qualifiers.keys():
                    if 'gene' in genes.qualifiers.keys():
                        gene_name=genes.qualifiers["gene"][0]
                    else:
                        gene_name=genes.qualifiers["locus_tag"][0]    
                    gene_locustag=genes.qualifiers["locus_tag"][0]#genes IDS
                    gene_extract=genes.extract(seq_object)
                    gene_extract.id=gene_name
                    gene_seq=gene_extract.seq #genes sequences
                    gene_length=len(gene_seq) #genes length
                    gene_gc=GC(gene_seq) # genes GC%
                    # counting genes nucleotides
                    Adinine=gene_seq.count("A")
                    Thymine=gene_seq.count("T")
                    Cytosine=gene_seq.count("C")
                    Guanine=gene_seq.count("G")
                    gene_mw=g_mw=mw(gene_seq,"DNA")/1000 #gene molecular weight
                    st_Tm_GC="%0.2f" % mt.Tm_GC(gene_seq)#Standard GC Melting Temperature
                    st_Tm_NN="%0.2f" % mt.Tm_NN(gene_seq)#Standard NN Melting Temperature
                    #data lists
                    gene_info=(gene_locustag,gene_name,gene_length,Adinine,Thymine,Cytosine,Guanine,gene_gc)
                    genes_mw=(gene_locustag,gene_name,gene_length,gene_mw)
                    genes_mt=(gene_locustag,gene_name,st_Tm_GC,st_Tm_NN)
                    #append lists
                    genes_locustages.append(gene_locustag)
                    genes_names.append(gene_name)
                    genes_info.append(gene_info)
                    genes_sequences.append(gene_extract)
                    genes_molecular_weight.append(genes_mw)
                    st_melting_temp.append(genes_mt)
                    #Prepare download all gene
                    g_seq=str(gene_seq)
                    n=60
                    splitted=[g_seq[i:i+n] for i in range(0,len(g_seq),n)]
                    all_g_seq_reslt='\n'.join(splitted)
                    g_desc=gene_extract.description
                    g_id_seq='>'+gene_name+' '+g_desc+'\n'+all_g_seq_reslt
                    all_g_fa.append(g_id_seq)
            #genomics dataframe
            genes_info_df=pd.DataFrame.from_records(genes_info,columns=["Locus Tag ID",'Gene ID',"Gene Length(bp)",'Adinine (bp)','Thymine (bp)','Cytosine (bp)','Guanine (bp)',"GC Content(%)",])
            genes_mw_df=pd.DataFrame.from_records(genes_molecular_weight,columns=['Locus Tag ID','Gene ID','Gene Length(bp)','Estimated MW (KDa)'])
            melting_temp_df=pd.DataFrame.from_records(st_melting_temp,columns=['Locus Tag ID','Gene ID','Standard Tm_GC (°C)','Standard Tm_NN (°C)'])
            #genes info
            st.write("**1- Genes Info**")
            st.markdown('')
            #filter genes
            filter_1=st.radio('Filter: ',["All Genes","Per Gene ID"],horizontal=True)
            if filter_1=="All Genes":
                st.dataframe(genes_info_df)
                st.markdown('')
                #download all genes_sequences
                allgenes_string='\n'.join(map(str,all_g_fa))
                st.write('Download all genes sequences:')
                st.download_button(label='Download Fasta',data=allgenes_string,file_name='all_genes.fasta',mime='text/plain')
            else:
                genes_filters=st.multiselect("Select Gene ID(s)",genes_names)
                filtred_genes_info=genes_info_df[genes_info_df["Gene ID"].isin(genes_filters)]
                st.dataframe(filtred_genes_info)
                # download selected genes
                st.markdown('')
                selected_genes_fasta=[]
                for record in genes_sequences:
                    gene_ID=record.id     
                    if gene_ID in genes_filters:
                        g_id=gene_ID
                        s_gene_seq=record.seq
                        s_g_seq=str(s_gene_seq)
                        n=60
                        s_splitted=[s_g_seq[i:i+n] for i in range(0,len(s_g_seq),n)]
                        s_g_seq_reslt='\n'.join(s_splitted)
                        s_g_description=record.description 
                        sg_id_seq='>'+g_id+' '+s_g_description+'\n'+s_g_seq_reslt
                        selected_genes_fasta.append(sg_id_seq)
                sg_string='\n'.join(map(str,selected_genes_fasta))
                st.write('Download selected gene(s):')
                st.download_button(label='Download Fasta',data=sg_string,file_name='selected_genes.fasta',mime='text/plain')
            st.markdown('')
            ##Physico-Chemical Properties
            st.write("**2- Physico-Chemical Properties**")
            st.markdown('')
            st.write("*2-1 Estimated Molecular Weight*")
            #filter genes
            filter_2=st.radio('Filter MW: ',['All Genes','Per Gene ID'],horizontal=True)
            if filter_2=='All Genes':
                st.dataframe(genes_mw_df)
            else:
                mw_filters=st.multiselect("Select Gene(s)",genes_names)
                filtred_mw=genes_mw_df[genes_mw_df["Gene ID"].isin(mw_filters)]
                st.dataframe(filtred_mw)
            ##Melting Temperature########
            st.markdown('')
            st.write('*2-2 Melting Temperature*')
            filter_3=st.radio('Filter MT: ',['All Genes','Per Gene ID'],horizontal=True)
            if filter_3=='All Genes':
                st.dataframe(melting_temp_df)
                st.write('**Tm_NN**: Melting.T° value based on nearest neighbor thermodynamics.')
                st.write('**Tm_GC**: Melting temperature value based on GC content.')
            
            else:
                mt_filters=st.multiselect("",genes_names)
                filtred_mt=melting_temp_df[melting_temp_df["Gene ID"].isin(mt_filters)]
                st.dataframe(filtred_mt)
                st.write('**Tm_NN**: Melting.T° value based on nearest neighbor thermodynamics.')
                st.write('**Tm_GC**: Melting temperature value based on GC content.')
            #########Proteomics#########
            st.markdown('***')
            st.header('Proteomics')
            #extract CDS
            proteins_features=[]
            for feature in seq_object.features:
                if feature.type=="CDS":
                    proteins_features.append(feature)
            #extract proteins
            #protein empty list
            p_products=[]
            genes_proteins=[]
            p_info_dataframe=[]
            ps_MW=[]
            p_AII=[]
            ps_MEC=[]
            ps_SSF=[]
            all_p_fa=[]
            for protein in proteins_features:
                if "translation" in protein.qualifiers.keys():
                    p_name=protein.qualifiers['product'][0]
                    p_id=protein.qualifiers['protein_id'][0]
                    protein_seq=protein.qualifiers["translation"][0]
                    p_length=len(protein_seq)
                    amino_acids=str(protein_seq)
                    #analysed Proteins
                    proteins_analysed=ProteinAnalysis(amino_acids)
                    #Mol Weight
                    p_MW="%0.2f" % proteins_analysed.molecular_weight()
                    # Iso Point,Arm, I.index and Gravy
                    p_IP="%0.2f" % proteins_analysed.isoelectric_point()
                    p_Arm="%0.2f" % proteins_analysed.aromaticity()
                    p_Idx="%0.2f" % proteins_analysed.instability_index()
                    p_GV="%0.2f" % proteins_analysed.gravy()
                    #molar extinction coefficient
                    p_MEC=proteins_analysed.molar_extinction_coefficient()
                    reduced_cysteines=p_MEC[0]
                    disulfid_bridges=p_MEC[1]
                    # Secondary Stracture Fractions
                    p_SSF=proteins_analysed.secondary_structure_fraction()
                    Helix="%0.2f" % p_SSF[0]
                    Turn="%0.2f" % p_SSF[1]
                    Sheet="%0.2f" % p_SSF[2]
                    #all protein seqs download preparation
                    n=60
                    splitted=[amino_acids[i:i+n] for i in range(0,len(amino_acids),n)]
                    all_p_seq_reslt='\n'.join(splitted)
                    allprotein_fasta='>'+p_id+' '+p_name+'\n'+all_p_seq_reslt
                    if "gene" in protein.qualifiers.keys():
                        p_gene=protein.qualifiers["gene"][0]
                    else:
                        p_gene=protein.qualifiers["locus_tag"][0]
                    #data
                    p_info=(p_name,p_id,p_gene,p_length)
                    ps_Mol_weight=(p_name,p_id,p_gene,p_MW)
                    AII=(p_name,p_id,p_gene,p_IP,p_Arm,p_Idx,p_GV)
                    MEC=(p_name,p_id,p_gene,reduced_cysteines,disulfid_bridges)
                    SSF=(p_name,p_id,p_gene,Helix,Turn,Sheet)
                    # append lists
                    p_products.append(p_name)
                    genes_proteins.append(p_gene)
                    p_info_dataframe.append(p_info)
                    ps_MW.append(ps_Mol_weight)
                    p_AII.append(AII)
                    ps_MEC.append(MEC)
                    ps_SSF.append(SSF)
                    all_p_fa.append(allprotein_fasta) #append download
            #proteins DataFrames
            p_info_df=pd.DataFrame(p_info_dataframe,columns=["Product","Protein ID","Gene ID","Protein Length(aa)"])
            ps_MW_df=pd.DataFrame(ps_MW,columns=["Product","Protein ID","Gene ID","Molecular Weight(Da)"])
            p_AII_df=pd.DataFrame(p_AII,columns=["Product","Protein ID","Gene ID","Isoelectric Point","Aromaticity","Instability Index","GRAVY"])
            ps_MEC_df=pd.DataFrame(ps_MEC,columns=["Product","Protein ID","Gene ID","MEC with reduced_cysteines","MEC with disulfid_bridges"])
            ps_SSF_df=pd.DataFrame(ps_SSF,columns=["Product","Protein ID","Gene ID","Helix","Turn","Sheet"])
            #proteins Info
            st.markdown('')
            st.subheader('I- Proteins Info')
            st.markdown('')
            st.write("**1- Proteins Data**")
            #filter
            p_info_filter=st.radio('',['All Proteins','Per Protein'],horizontal=True)
            if p_info_filter=='All Proteins':
                st.dataframe(p_info_df)
                allproteins_string='\n'.join(map(str,all_p_fa))
                st.write('Download all proteins sequences:')
                st.download_button(label='Download Fasta',data=allproteins_string,file_name='all_proteins.fasta',mime='text/plain')
            else:
                p_f_p=st.multiselect('Select Product(s)',p_products)
                p_p_filtered=p_info_df[p_info_df["Product"].isin(p_f_p)]
                st.dataframe(p_p_filtered)
                #download selected protein sequences
                s_p_f=[]
                for protein in proteins_features:
                    if "translation" in protein.qualifiers.keys():
                        p_name=protein.qualifiers['product'][0]
                        if p_name in p_f_p:
                            p_id=protein.qualifiers['protein_id'][0]
                            protein_seq=protein.qualifiers["translation"][0]
                            amino_acids=str(protein_seq)
                            n=60
                            s_splitted=[amino_acids[i:i+n] for i in range(0,len(amino_acids),n)]
                            s_p_seq_reslt='\n'.join(s_splitted)
                            sp_id_seq='>'+p_id+' '+p_name+'\n'+s_p_seq_reslt
                            s_p_f.append(sp_id_seq)
                sp_string='\n'.join(map(str,s_p_f))
                st.markdown('')
                st.write('Download selected protein(s):')
                st.download_button(label='Download Fasta',data=sp_string,file_name='selected_protein.fasta',mime='text/plain')
            st.markdown('')
            st.subheader('II- Physico-Chemical Properties')
            st.markdown('')
            st.write("**1- Estimated Molecular Weight**")
            st.markdown("")
            p_mw_filter=st.radio("Filter Protein (MW)",['All Proteins','Per Protein'],horizontal=True)
            if p_mw_filter=="All Proteins":
                st.dataframe(ps_MW_df)
            else:
                pmw_f_p=st.multiselect("Select Product",p_products)
                pmw_p_filtered=ps_MW_df[ps_MW_df["Product"].isin(pmw_f_p)]
                st.dataframe(pmw_p_filtered)
            st.markdown('')
            st.write("**2- Aromaticity (Lobry,1994), Isoelectric Point, Instability Index(Guruprasad et al.1990) & GRAVY (Kyte and Doolittle,1982)**")
            st.markdown('')
            p_aii_filter=st.radio("Filter(AIIG)",['All Proteins','Per Protein'],horizontal=True)
            if p_aii_filter=="All Proteins":
                st.dataframe(p_AII_df)
            else:
                paii_f_p=st.multiselect("Select",p_products)
                paii_p_filtered=p_AII_df[p_AII_df["Product"].isin(paii_f_p)]
                st.dataframe(paii_p_filtered)
            # Estimate the molar extinction coefficient
            st.markdown('')
            st.write("**3- Estimated Molar Extinction Coefficient**")
            st.markdown("")
            p_mec_filter=st.radio("Filter(MEC)",['All Proteins','Per Protein'],horizontal=True)
            if p_mec_filter=="All Proteins":
                st.dataframe(ps_MEC_df)
            else:
                pmec_f_p=st.multiselect("",p_products)
                pmec_p_filtered=ps_MEC_df[ps_MEC_df["Product"].isin(pmec_f_p)]
                st.dataframe(pmec_p_filtered)
            # Estimate the Secondary Stracture Fractions
            st.markdown('')
            st.write("**4- Secondary Structure Fractions**")
            st.markdown('')
            p_ssf_filter=st.radio("Filter(SSF)",['All Proteins','Per Protein'],horizontal=True)
            if p_ssf_filter=="All Proteins":
                st.dataframe(ps_SSF_df)
            else:
                pssf_f_p=st.multiselect("Select Protein",p_products)
                pssf_p_filtered=ps_SSF_df[ps_SSF_df["Product"].isin(pssf_f_p)]
                st.dataframe(pssf_p_filtered)
            
        except IncompleteRead:
            st.write('Check Your Connexion And Try Again')
    #Anammox Genra selection
    st.sidebar.markdown('***')
    select_genra=st.sidebar.selectbox('Select Genra',options=['Candidatus Kuenenia','Candidatus Jettenia','Candidatus Brocadia','Candidatus Scalindua', 'Ca. Anammoxomicrobium','Candidatus Anammoxoglobus'])
    
    if select_genra=='Candidatus Kuenenia':
        select_species=st.selectbox('Select Species',['Ca. Kuenenia Stuttgartiensis'])
        if select_species=='Ca. Kuenenia Stuttgartiensis':
            select_strain=st.selectbox('Select strain/Isolate',['Ca. Kuenenia stuttgartiensis MBR-1','Ca. Kuenenia stuttgartiensis CSTR1','Ca. Kuenenia stuttgartiensis Kaust'])
            if select_strain=='Ca. Kuenenia stuttgartiensis MBR-1':
                ANID_LTGN('LT934425')
            elif select_strain=='Ca. Kuenenia stuttgartiensis CSTR1':
                ANID_LTGN('CP049055')
            else:
                contig_list=['CT{}'.format(i) for i in range(573071,573076)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
    #Candidatus Jettenia#############            
    elif select_genra=='Candidatus Jettenia':
        select_species=st.selectbox('Select Species',['Ca. Jettenia Caeni','Ca. Jettenia Ecosi'])
        if select_species=='Ca. Jettenia Caeni':
            select_strain=st.selectbox('Select strain/Isolate',['Jettenia Caeni KSU1','Jettenia Caeni JETCAE04'])
            if select_strain=='Jettenia Caeni KSU1':
                contig_list=['BAFH0{}'.format(i) for i in range(1000001,1000005)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
            else:
                contig_list=['BQMT0{}'.format(i) for i in range(1000001,1000096)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
        else:
            select_strain=st.selectbox('Select strain/Isolate',['Jettenia Ecosi isolate J2'])
            if select_strain=='Jettenia Ecosi isolate J2':
                contig_list=['SULG0{}'.format(i) for i in range(1000001,1000224)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)           
    #Candidatus Brocadia#############
    elif select_genra=='Candidatus Brocadia':
        select_species=st.selectbox('Select Species',['Ca. Brocadia Fulgida','Ca. Brocadia Sinica','Ca. Brocadia Sapporoensis','Ca. Brocadia Caroliniensis'])
        if select_species=='Ca. Brocadia Fulgida':
            select_strain=st.selectbox('Select strain/Isolate',['Brocadia fulgida isolate RU1'])
            if select_strain=='Brocadia fulgida isolate RU1':
                contig_list=['LAQJ0{}'.format(i) for i in range(1000001,1000316)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
        elif select_species=='Ca. Brocadia Sinica':
            select_strain=st.selectbox('Select strain/Isolate',['Brocadia sinica isolate OLB1'])
            if select_strain=='Brocadia sinica isolate OLB1':
                contig_list=['JZEK0{}'.format(i) for i in range(1000001,1000087)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
        elif select_species=='Ca. Brocadia Sapporoensis':
            select_strain=st.selectbox('Select strain/Isolate',['Brocadia sapporoensis strain 40','Brocadia sapporoensis HBSAPP01'])
            if select_strain=='Brocadia sapporoensis strain 40':
                contig_list=['MJUW0{}'.format(i) for i in range(2000001,2000124)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
            else:
                contig_list=['BQMM0{}'.format(i) for i in range(1000001,1000140)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
        else:
            select_strain=st.selectbox('Select strain/Isolate',['Brocadia Caroliniensis Isolate 26THWARD'])
            if select_strain=='Brocadia Caroliniensis Isolate 26THWARD':
                contig_list=['AYTS0{}'.format(i) for i in range(1000001,1000210)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)   
    #Candidatus Scalindua############    
    elif select_genra=='Candidatus Scalindua':
        select_species=st.selectbox('Select Species',['Ca. Scalindua Japonica','Ca. Scalindua Arabica','Ca. Scalindua Brodae','Ca. Scalindua Sp'])
        if select_species=='Ca. Scalindua Japonica':
            select_strain=st.selectbox('Select strain/Isolate',['Scalindua japonica strain husup-a2'])
            if select_strain=='Scalindua japonica strain husup-a2':
                contig_list=['BAOS0{}'.format(i) for i in range(1000001,1000048)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
        elif select_species=='Ca. Scalindua Arabica':
            select_strain=st.selectbox('Select strain/Isolate',['Scalindua arabica isolate SuakinDeep_MAG55_1'])
            if select_strain=='Scalindua arabica isolate SuakinDeep_MAG55_1':
                contig_list=['JAANXD0{}'.format(i) for i in range(10000002,10000103)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
        elif select_species=='Ca. Scalindua Brodae':
            select_strain=st.selectbox('Select strain/Isolate',['Scalindua brodae isolate RU1'])
            if select_strain=='Scalindua brodae isolate RU1':
                contig_list=['JRYO0{}'.format(i) for i in range(1000001,1000283)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
        else:
            select_strain=st.selectbox('Select strain/Isolate',['Scalindua sp. SCAELEC01','Scalindua sp. AMX11'])
            if select_strain=='Scalindua sp. SCAELEC01':
                contig_list=['SHMT0{}'.format(i) for i in range(1000001,1000105)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
            else:
                contig_list=['RBMW0{}'.format(i) for i in range(1000001,1000122)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
        
    #Ca. Anammoxomicrobium############
    elif select_genra=='Ca. Anammoxomicrobium':
        select_species=st.selectbox('Select Species',['Ca. Anammoximicrobium sp'])
        if select_species=='Ca. Anammoximicrobium sp':
            select_strain=st.selectbox('Select strain/Isolate',['Anammoximicrobium sp. isolate AS06rmzACSIP_251 588_AS06'])
            if select_strain=='Anammoximicrobium sp. isolate AS06rmzACSIP_251 588_AS06':
                contig_list=['JAAYYU0{}'.format(i) for i in range(10000001,10000418)]
                contig_menu=st.selectbox('Select Contig',contig_list)
                ANID_LTGN(contig_menu)
    
    #Ca. Anammoxoglobus Propionicus####
    else:
        select_species=st.selectbox('Select Species',['Ca. Anammoxoglobus Propionicus'])
        st.markdown('')
        st.write('Complete genome is not available')
        
#############################################################################################################################################References
############################################################################################################################################
else:
    st.sidebar.markdown('***')
    #about
    st.sidebar.title('About')
    st.sidebar.info('''
    This app is made using :\n 
    Streamlit - Python - Biopython - Pandas- Py3dmol - Matplotlib - Leafmap\n\n
    It is maintained by [Yassine OUCHEN](https://github.com/ouchen-bio/Anammox-Bacteria-Database)
    Email: yassine.ouchen2021@gmail.com
    LinkedIn: [Click Here](https://www.linkedin.com/in/yassine-ouchen-be)
    ''')
    st.subheader('References')
    st.markdown('***')
    st.write("**Akram M., Dietl A., Mersdorf U., et al. (2015)**. *A 192-heme electon transfer network in the hydrazine dehydrogenase complex. Sci Adv **5**: eaav4310.*")
    st.markdown("")
    st.write("**Ali M., Oshiki M., Awata T., Isobe K.,Kimura, Yoshikawa H., Hira D., Kindaichi T., Satoh H., Fujii T., Okabe S. (2015)**. *Physiological characterization of anaerobic ammonium oxidizing bacterium ‘Candidatus/Jettenia caeni’. Environ. Microbiol. **17**:  2172–2189.*")
    st.markdown("")
    st.write('**Arrojo B., Mosquera-Corral A. (2006)**. *Effects of mechanical stress on anammox granules in a sequencing batch reactor (SBR). J. Biotechnol. **123**:  453–463.*')
    st.markdown("")
    st.write('**Awata T., Oshiki M., Kindaichi T., Ozaki N., Ohashi A., and Okabe S. (2013)**. *Physiological characterization of an anaerobic ammonium-oxidizing bacterium belonging to the ‘Candidatus Scalindua’ group. Appl Environ Microbiol **79**: 4145–4148*.')
    st.markdown("")
    st.write("**Brandsma J., van de Vossenberg J., RisgaardPetersen N., Schmid M.C., Engstro¨m, P., Eurenius K., Hulth S., Jaeschke A., Abbas B., Hopmans E.C., Strous M., Schouten S., Jetten M.S.M., Damste´ J.S.S (2011)**. *A multi-proxy study of anaerobic ammonium oxidation in marine sediments of the Gullmar Fjord, Sweden, Environ. Microbiol. Rep. **3**: 360–366.*")
    st.markdown("")
    st.write('**Botchkova E. A., Litti Yu. V., Grouzdev D., Novikov A., Bochkareva E.S., Beskorovayny A.V., Kuznetsov B.B., Nozhevnikova A.N. (2018)**. *Description of “Candidatus Jettenia ecosi” sp. nov., a New Species of Anammox Bacteria. Microbiology. **87**: 766-776.*')
    st.markdown('')
    st.write('**Casagrande C.G., Kunz A., de Pra´M.C., Bressan C.R., Soares H.M. (2013)**. *High nitrogen removal rate using anammox process at short hydraulic retention time, Water Sci. Technol. **67**: 968–975.*')
    st.markdown("")
    st.write('**Carvajal-Arroyo J.M., Sun W., Sierra-Alvarez R., and Field J.A. (2013)**. *Inhibition of anaerobic ammonium oxidizing (anammox) enrichment cultures by substrates, metabolites and common wastewater constituents. Chemosphere **91**: 22–27.*')
    st.markdown("")
    st.write('**Chen J., Zheng P., Yu Y., Tang C., Mahmood Q. (2010)**. *Promoting sludge quantity and activity results in high loading rates in anammox UBF, Bioresour. Technol. **101**:  2700–2705.*')
    st.markdown("")
    st.write("**Chicano T.M., Dietrich L., de Almeida N.M, Akram M., Hartmann E., Leidreiter F., Leopoldus D., Mueller M., Sánchez R., Nuijten G.H.L., Reimann J., Seifert K.A, Schlichting I., van Niftrik L., Jetten M.S.M., Dietl A., Kartal B., Parey K., Barends T.R.M. (2021)**. *Structural and functional characterization of the intracellular filament-forming nitrite oxidoreductase multiprotein complex. Nature Microbiology VOL **6**: 1129–1139 .*")
    st.markdown("")
    st.write('**Damste J.S.S., Strous M., Rijpstra W.I.C., Hopmans E.C., Geenevasen J.A.J., van Duin A.C.T., van Niftrik L., Jetten M.S.M. (2002)**. *Linearly concatenated cyclobutane lipids form a dense bacterial membrane. Nature.**419**: 708-712.*')
    st.markdown("")
    st.write('**Dapena-Mora A., Campos J.L., Mosquera-Corral A., Jetten M.S.M., Mendez R. (2004)**. *Stability of the anammox process in a gas-lift reactor and a SBR. J. Biotechnol. **110**:  159–170.*')
    st.markdown("")
    st.write("**Egli K., Fanger U., Alvarez P.J.J., Siegrist H., van der Meer J.R., Zehnder A.J.B. (2001)**. *Enrichment and characterization of an anammox bacterium from a rotating biological contactor treating ammonium-rich leachate, Arch. Microbiol. **175**: 198–207.*")
    st.markdown("")
    st.write("**Fuerst JA. (2005)**. *Intracellular compartmentation in planctomycetes. Annu Rev Microbiol **59**: 299–328.*")
    st.markdown("")
    st.write("**Fuerst JA. & Sagulenko E. (2011)**. *Beyond the bacterium: planctomycetes challenge our concepts of microbial structure and function. Nat Rev Microbiol **9**: 403–413.*")
    st.markdown("")
    st.write('**Gori F., Tringe S.G., Kartal B., Machiori E., and Jetten M.S.M. (2011)**. *The metagenomic basis of anammox metabolism in Candidatus ‘Brocadia fulgida. Biochem Soc Trans **39**: 1799–1804*.')
    st.markdown("")
    st.write("**Guruprasad K., Reddy B.V., Pandit M.W. (1990)**. *Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence. Protein Eng. **4 (2)**: 155–61.*")
    st.markdown("")
    st.write('**Hooper A.B., Vannelli T., Bergmann D.J., Arciero D.M. (1997)**. *Enzymology of the oxidation of ammonia to nitrite by bacteria. Antonie Van Leeuwenhoek **71**: 59-67.*')
    st.markdown("")
    st.write('**Hopmans E.C., Kienhuis M.V.M., Rattray J.E., Jaeschke A., Schouten S., and Sinninghe Damsté J.S. (2006)**. *Improved analysis of ladderane lipids in biomass and sediments using high-performance liquid chromatography/atmospheric pressure chemical ionization tandem mass spectrometry. Rapid Communications in Mass Spectrometry **20**: 2099-2103.*')
    st.markdown('')
    st.write('**Hu Z., van Alen T., Jetten M.S., and Kartal B. (2013b)**. *Effects of lysozyme and penicillin on the growth and activity of anaerobic ammonium-oxidizing planctomycetes. Appl Environ Microbiol **79**: 7763–7769.*')
    st.markdown('')
    st.write("**Jetten M.S.M., Logemann S., Muyzer G., Robertson L.A., de Vries S., van Loosdrecht M.C. and Kuenen J.G. (1997)**. Novel principles in the microbial conversion of nitrogen compounds. Antonie Van Leeuwenhoek **71**: 75–93.")
    st.markdown("")
    st.write('**Jetten M.S.M., Niftrik L., Strous M., Kartal B., Keltjens J.T., Op den Camp H.J. (2009)**. *Biochemistry and molecular biology of anammox bacteria. Critical Reviews in Biochemistry and Molecular Biology **44**, 65e84*.')
    st.markdown("")
    st.write('**Jiang T., Zhang H., Qiang H., Yang F., Xu X., Du H. (2013)**. *Start-up of the anammox process and membrane fouling analysis in a novel rotating membrane bioreactor, Desalination **311**: 46–53*.')
    st.markdown("")
    st.write('**Jin R.C., Zheng P., Hu A.H., Mahmood Q., Hu B.L., Jilani G. (2008)**. *Performance comparison of two anammox reactors: SBR and UBF. Chem. Eng. J. **138**:  224–230.*')
    st.markdown("")
    st.write('**Jin R.C., Ma C., Mahmood Q., Yang G.F., Zheng P. (2011)**. *Anammox in a UASB reactor treating saline wastewater, Process Saf. Environ. Prot. **89**: 342–348.*')
    st.markdown("")
    st.write("**Kartal B., Rattray J., van Niftrik L. et al. (2007a)**. *Candidatus “Anammoxoglobus propionicus” a new propionate oxidizing species of anaerobic ammonium oxidizing bacteria. Syst Appl Microbiol **30**: 39–49.*")
    st.markdown("")
    st.write("**Kartal B., van Niftrik L., Rattray J., van de Vossenberg J.L.C.M., Schmid M.C., Damste J.S., Jetten M.S.M., Strous M. (2008)**. *Candidatus ‘Brocadia fulgida’: an autofluorescent anaerobic ammonium oxidizing bacterium, FEMS Microbiol. Ecol. **63**: 46–55.*")
    st.markdown("")
    st.write('**Kartal B., van Niftrik L., Keltjens J.T., Op den Camp H.J.M., and Jetten M.S.M. (2012)**. *Anammox–growth physiology, cell biology, and metabolism. Adv Microb Physiol **60**: 211–262.*')
    st.markdown("")
    st.write('**Kartal B., de Almeida N.M., Maalcke W.J., Op den Camp H.J.M., Jetten M.S.M., and Keltjens J.T. (2013)**. *How to make a living from anaerobic ammonium oxidation. FEMS Microbiol Rev **37**: 428–461.*')
    st.markdown("")
    st.write("**Khramenkov S.V., Kozlov M.N., Kevbrina M.V., Dorofeev A.G., Kazakova E.A., Grachev V.A., Kuznetsov B.B., Polyakov D.Y., Nikolaev Y.A. (2013)**. *A novel bacterium carrying out anaerobic ammonium oxidation in a reactor for biological treatment of the filtrate of wastewater fermented sludge, Microbiology **82**: 628–636.*")
    st.markdown("")
    st.write('**Kindaichi T., Awata T., Tanabe K., Ozaki N., Ohashi A. (2011)**. *Enrichment of marine anammox bacteria in Hiroshima Bay sediments, Water Sci. Technol. **63**: 965–970.*')
    st.markdown("")
    st.write('**Kuypers M.M.M., Sliekers A.O., Lavik G., Schmid M., Jørgensen B.B., Kuenen J.G. et al; (2003)**. *Anaerobic ammonium oxidation by anammox bacteria in the Black Sea. Nature **422**: 608-611.*')
    st.markdown('')
    st.write('**Kuypers M.M.M., Lavik G., Woebken D., Schmid M., Fuchs B.M., Amann R. et al. (2005)**. *Massive nitrogen loss from the Benguela upwelling system through anaerobic ammonium oxidation. Proc Natl Acd Sci USA **102**: 6478-6483.*')
    st.markdown('')
    st.write("**Kyte J. and Doolittle R. (1982).** *A simple method for displaying the hydropathic character of a protein. J. Mol. Biol. **157**: 105-132.*")
    st.markdown("")
    st.write('**Lam P., Lavik G., Jensen M.M., van de Vossenberg J., Schmid M., Woebken D., Gutierrez D., Amann R., Jetten M.S.M., Kuypers M.M. (2009)**. *Revising the nitrogen cycle in the Peruvian oxygen minimum zone. Proc Natl Acad Sci USA **106**:4752–4757.*')
    st.markdown('')
    st.write("**Lindsay M.R., Webb R.I., Strous M., Jetten M.S.M., Butler M.K., Forde R.J. and Fuerst J.A. (2001)**. *Cell compartmentalisation in planctomycetes: novel types of structural organisation for the bacterial cell. Arch. Microbiol. **175**: 413–429.*")
    st.markdown("")
    st.write('**Li M., Ford T., Li XY., Gu J-D. (2011b)**. *Cytochrome cd-1-containing nitrite encoding gene nirS as a new functional biomarker for detection of anaerobic ammonium oxidizing (Anammox) bacteria. Envion Sci Technol. doi:10.1021/es103826w.*')
    st.markdown('')
    st.write('**Li H., Zhou S., Ma W., Huang G., Xu B. (2012)**. *Fast start-up of anammox reactor: operational strategy and some characteristics as indicators of reactor performance, Desalination **286**: 436–441.*')
    st.markdown("")
    st.write("**Li Z., Xu X., Yang F., Zhang S. (2015)**. *Sustainable operation of submerged Anammox membrane bioreactor with recycling biogas sparging for alleviating membrane fouling, Chemosphere. **140**: 106–113*.")
    st.markdown("")
    st.write("**Lobry J.R., and Gautier C. (1994)**. *Hydrophobicity, expressivity and aromaticity are the major trends of amino acid usage in 999 Escherichia coli chromosome encoded genes. Nucleic Acids Research **22**: 3174-3180.*")
    st.markdown("")
    st.write("**Mulder A., van de Graaf A. A., Robertson L. A. & Kuenen J. G. (1995)**. *Anaerobic ammonium oxidation discovered in a denitrifying fluidized bed reactor. FEMS Microbiol Ecol **16**: 177–184.*")
    st.markdown("")
    st.write('**Murray R.G.E. and Stackebrandt E. (1995)**. *Taxonomic note: implementation of the provisional status Candidatus for incompletely described procaryotes, Int. J. Syst. Bacteriol. vol. **45**, pp. 186–187.*')
    st.markdown('')
    st.write('**Narita Y., Zhang L., Ali M., Kimura Z.I. (2017)**. *Enrichment and physiological characterization of an anaerobic ammonium-oxidizing bacterium ‘Candidatus Brocadia sapporoensis’. syst appl microbiol **40**: 448-457*')
    st.markdown("")
    st.write('**Ni S.Q., Gao B.Y., Wang C.C., Lin J.G., Sung S. (2011)**. *Fast start-up, performance and microbial community in a pilot-scale anammox reactor seeded with exotic mature granules, Bioresour. Technol. **102**: 2448–2454.*')
    st.markdown('')
    st.write('**Nicholas Rego and David Koes (2015)**. *3Dmol.js: molecular visualization with WebGL. Bioinformatics. **31**: 1322-1324.*')
    st.markdown('')
    st.write('**Oshiki M., Shimokawa M., Fujii N., Satoh H., and Okabe S. (2011)**. *Physiological characteristics of the anaerobic ammonium-oxidizing bacterium ‘Candidatus Brocadia sinica. Microbiology **157**: 1706–1713*.')
    st.markdown("")
    st.write('**Oshiki M., Awata T., Kindaichi T., Satoh H., and Okabe S. (2013a)**. *Cultivation of planktonic anaerobic ammonium oxidation (anammox) bacteria by using membrane bioreactor. Microbes Environ **28**: 436–443*.')
    st.markdown("")
    st.write('**Oshiki M., Ishii S., Yoshida K., Fujii N., Ishiguro M., Satoh H., et al. (2013b)**. *Nitrate-dependent ferrous iron oxidation by anaerobic ammonium oxidation (anammox) bacteria. Appl Environ Microbiol **79**: 4087–4093*.')
    st.markdown("")
    st.write('**Puyol D., Carvajal-Arroyo J.M., Garcia B., Sierra-Alvarez R., and Field J.A. (2013)**. *Kinetic characterization of Brocadia spp.-dominated anammox cultures. Bioresour Technol **1391**: 94–100.*')
    st.markdown("")
    st.write("**Quan Z.X, Rhee SK., Zuo JE., Yang Y., Bae JW., Park JR., Lee ST. & Park YH. (2008)**. *Diversity of ammonium-oxidizing bacteria in a granular sludge anaerobic ammoniumoxidizing (anammox) reactor. Environ Microbiol **10**: 3130–3139.*")
    st.markdown("")
    st.write('**Rattray J., van Vossenberg J., Hopmans E., Kartal B., van Niftrik L., Rijpstra W., Strous M., Jetten M.S.M., Schouten S., Damste J.S.S. (2008)**. *Ladderane lipid distribution in four genra of anammox bacteria. Arch of Microbiology **190**: 51-66*.')
    st.markdown("")
    st.write("**Schmid M., Walsh K., Webb R. et al. (2003)**. *Candidatus “Scalindua brodae”, sp. nov., Candidatus “Scalindua wagneri”, sp. nov., two new species of anaerobic ammonium oxidizing bacteria. Syst Appl Microbiol **26**: 529–538.*")
    st.markdown("")
    st.write('**Schmid M.C, Hooper A.B, Klotz M.G., Woebken D., Lam P., Kuypers M.M., Pommerening-Roeser A., Op den Camp HJ, Jetten M.S.M (2008)**. *Environmental detection of octahaem cytochrome c hydroxylamine/hydrazine oxidoreductase genes of aerobic and anaerobic ammonium-oxidizing bacteria. Environ Microbiol **10**:3140–3149.*')
    st.write('**Sneath P.H.A. (1990)**. *International Code of Nomenclature of Bacteria. Revision, Washington: ASM, 1990.*')
    st.markdown("")
    st.write('**Strous M., Heijnen J.J., Kuenen J.G., and Jetten M.S.M. (1998)**. *The sequencing batch reactor as a powerful tool for the study of slowly growing anaerobic ammonium-oxidizing microorganisms. Appl Microbiol Biotechnol **50**: 589–596.*')
    st.markdown("")
    st.write("**Strous M., Fuerst J.A., Kramer E.H., Logemann S., Muyzer G., van de Pas-Schoonen K.T., Webb R., Kuenen J.G. & Jetten M.S.M. (1999a)**. *Missing lithotroph identified as new planctomycete. Nature **400**: 446–449.*")
    st.markdown("")
    st.write('**Strous M., Gijs Kuenen J., and Jetten M.S.M. (1999b)**. *Key physiology of anaerobic ammonium oxidation. Appl Environ Microbiol **65**: 3248–3250.*')
    st.markdown("")
    st.write('**Strous M., Pelletier E., Mangenot S., Rattei T., Lehner A., Taylor M., et al. (2006)**. *Deciphering the evolution and metabolism of an anammox bacterium from a community genome. Nature **440**: 790–794.*')
    st.markdown("")
    st.write('**Trigo C., Campos J.L., Garrido J.M., Me´ndez R. (2006)**. *Startup of the anammox process in a membrane bioreactor. J. Biotechnol. **126**: 475–487.*')
    st.markdown("")
    st.write('**Van de Graaf A.A, de Bruijn P, Robertson L.A, Jetten M.S.M, Kuenen J.G. (1996)**. *Autotrophic growth of anaerobic ammonium-oxidizing micro-organisms in a fluidized bed reactor, Microbiology **142**: 2187–2196.*')
    st.markdown("")
    st.write('**Van de Vossenberg J., Rattray J.E., Geerts W., Kartal B., van Niftrik L., van Donselaar E.G., Damste´ J.S.S., Strous M., Jetten M.S.M. (2008)**. *Enrichment and characterization of marine anammox bacteria associated with global nitrogen gas production, Environ. Microbiol. **10**: 3120–3129.*')
    st.markdown("")
    st.write("**Van de Vossenberg J., Woebken D., Maalcke W.J., Wessels H.J., Dutilh B.E., Kartal B., Janssen-Megens E.M., Roeselers G., Yan  J., Speth D., Gloerich J., Geerts W., van der Biezen E., Pluk W., Francoijs K.J., Russ L., Lam P., Malfatti S.A., Tringe S.G., Haaijer S.C.M., Op den Camp H.J.M., Stunnenberg H.G., Amann R., Kuypers M.M.M., Jetten, M.S.M. (2013)**. *The metagenome of the marine anammox bacterium ‘Candidatus Scalindua profunda’ illustrates the versatility of this globally important nitrogen cycle bacterium, Environ. Microbiol. **15**: 1275–1289.*")
    st.markdown("")
    st.write('**Van der Star W.R.L., Abma W.R., Blommers D., Mulder J.W., Tokutomi T., Strous M., Picioreanu C., van Loosdrecht M.C. (2007)**. *Startup of reactors for anoxic ammonium oxidation: Experiences from the first full-scale anammox reactor in Rotterdam, Water Res. **41**: 4149–4163.*')
    st.markdown('')
    st.write("**Van der Star W.R.L., Miclea A.I., van Dongen U.G.J.M., Muyzer G., Picioreanu C., van Loosdrecht M.C.M. (2008)**. *The membrane bioreactor: a novel tool to grow anammox bacteria as free cells, Biotechnol. Bioeng. **101**: 286–294.*")
    st.markdown("")
    st.write("**Van Niftrik L., Geerts W.J., van Donselaar E.G., Humbel B.M., Yakushevska A., Verkleij A.J., Jetten M.S.M. and Strous M. (2008b)**. *Combined structural and chemical analysis of the anammoxosome: a membrane-bounded intracytoplasmic compartment in anammox bacteria. J. Struct. Biol. **161**: 401–410.*")
    st.markdown("")
    st.write("**Van Niftrik L. & Jetten M.S.M. (2012)**. *Anaerobic ammoniumoxidizing bacteria: unique microorganisms with exceptional properties. Microbiol Mol Biol Rev **76**: 585–596.*")
    st.markdown("")
    st.write("**Vera L., Gonza´lez E., Dı´az O., Delgado S. (2014)**. *Application of a backwashing strategy based on transmembrane pressure set-point in a tertiary submerged membrane bioreactor, J. Membr. Sci. **470**: 504–512.*")
    st.markdown("")
    st.write('**Wang T., Zhang H., Gao D., Yang F., Yang S., Jiang T., Zhang G. (2011)**. *Enrichment of anammox bacteria in seed sludges from different wastewater treating processes and start-up of anammox process, Desalination **271**: 193–198.*')
    st.markdown("")
    st.write("**Woebken D., Lam P., Kuypers M.M.M., Naqvi S.W.A., Kartal B., Strous M., Jetten M.S.M., Fuchs B.M., Amann R. (2008)**. *A microdiversity study of anammox bacteria reveals a novel Candidatus Scalindua phylotype in marine oxygen minimum zones, Environ. Microbiol. **10**: 3106–3119.*")
    st.markdown("")
    st.write('**Zhao R., Zhang H., Li Y., Jiang T., and Yang F. (2014)**. *Research of iron reduction and the iron reductase localization of anammox bacteria. Curr Microbiol **69**: 880–887.*')
