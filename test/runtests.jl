# this tests the julia interface to MAGEMin
using Test

cd("../")       # change to main directory

using MAGEMin_C

# Initialize database 
gv, DB = init_MAGEMin();

ble struct outP{ _T  } 
    P           ::  _T
    T           ::  _T 
    test        ::  Int64

    G           ::  _T
    ph          ::  Vector{String}
    ph_frac     ::  Vector{Float64}
end

# load reference for built-in tests
test0list       = Vector{outP}(undef, 42)
test0list       = outP[outP{Float64}(0.0, 800.0, 0, -810.5609033097174, ["ol", "cpx", "pl4T", "spn", "opx", "ru"], [0.667132274237546, 0.061999513468373686, 0.09433281071880929, 0.0021400860858959087, 0.1733594120033265, 0.0010359034860486904]), outP{Float64}(0.0, 1000.0, 0, -832.7707145063879, ["pl4T", "ol", "cpx", "spn", "opx", "spn"], [0.08839650154187155, 0.6611832051681615, 0.06169947650767187, 0.00026165843637477253, 0.1862452439333256, 0.0022139144125945852]), outP{Float64}(0.0, 1200.0, 0, -857.1069503923619, ["ol", "spn", "pl4T", "liq", "opx"], [0.7328729540247891, 0.00266101690764698, 0.011617667900985793, 0.21913853028003222, 0.03370066183597559]), outP{Float64}(0.0, 1400.0, 0, -883.9061639270158, ["spn", "ol", "liq"], [0.001693994239492248, 0.6964810512291004, 0.30182353927287503]), outP{Float64}(0.0, 1600.0, 0, -912.7435789030362, ["liq", "ol"], [0.43604649626734904, 0.563923059774354]), outP{Float64}(0.0, 1800.0, 0, -944.1811750598771, ["liq", "ol"], [0.9652440598206766, 0.034769170573530274]), outP{Float64}(0.0, 2000.0, 0, -978.6274324990358, ["liq"], [0.9999918060095452]), outP{Float64}(10.0, 800.0, 0, -794.613615369039, ["cpx", "opx", "ol", "spn"], [0.14140981969314273, 0.24252951850458013, 0.5873008276423136, 0.028759834159963336]), outP{Float64}(10.0, 1000.0, 0, -816.6838029287281, ["cpx", "opx", "spn", "ol"], [0.14694499013890033, 0.23699434805882247, 0.018303192176422185, 0.5977574696258551]), outP{Float64}(10.0, 1200.0, 0, -840.8076914926885, ["pl4T", "opx", "cpx", "ol", "spn", "liq"], [0.004662515046105702, 0.22245934622957136, 0.14864908360388715, 0.6127574541814642, 0.005801384045096177, 0.005669741069513012]), outP{Float64}(10.0, 1400.0, 0, -867.0488109158783, ["liq", "ol", "opx"], [0.19851781352106188, 0.6179054860717401, 0.1835798029379549]), outP{Float64}(10.0, 1600.0, 0, -895.561452915326, ["liq", "ol"], [0.4041308154005303, 0.5958689694876327]), outP{Float64}(10.0, 1800.0, 0, -926.4501402948199, ["liq", "ol"], [0.7464115581948364, 0.2535566570795361]), outP{Float64}(10.0, 2000.0, 0, -960.4694463085308, ["liq"], [0.9999969639994634]), outP{Float64}(20.0, 800.0, 0, -779.0911784097865, ["cpx", "g", "opx", "ol", "ru"], [0.11711774895821266, 0.12542453266284959, 0.14042422822002906, 0.6167416416519192, 0.000291848506989434]), outP{Float64}(20.0, 1000.0, 0, -801.016559084951, ["ol", "g", "opx", "cpx"], [0.6160606618022769, 0.11016691031456696, 0.15230898301409146, 0.12146344486906444]), outP{Float64}(20.0, 1200.0, 0, -824.9827248112628, ["cpx", "g", "opx", "ol"], [0.14058576833603953, 0.07457307374441834, 0.16878049611726487, 0.6160606618022771]), outP{Float64}(20.0, 1400.0, 0, -850.7858043518563, ["cpx", "liq", "ol", "opx"], [0.16324525231611625, 0.02983696105575796, 0.611553431422399, 0.19536448697422046]), outP{Float64}(20.0, 1600.0, 0, -878.8084720149644, ["liq", "ol", "opx"], [0.28968930948024535, 0.5520580980848422, 0.15826748027634258]), outP{Float64}(20.0, 1800.0, 0, -909.2903837645354, ["liq", "ol"], [0.6560433554673131, 0.3439261487009854]), outP{Float64}(20.0, 2000.0, 0, -942.8664583203127, ["liq"], [0.999989920442628]), outP{Float64}(30.0, 800.0, 0, -763.7568467920446, ["cpx", "g", "opx", "ol", "ru"], [0.11373256142153335, 0.14136671670432938, 0.12873766671952028, 0.616132337148915, 3.071800570207355e-5]), outP{Float64}(30.0, 1000.0, 0, -785.5631318330751, ["g", "opx", "cpx", "ol"], [0.13634836245717719, 0.13227049746136688, 0.11532047827917899, 0.6160606618022769]), outP{Float64}(30.0, 1200.0, 0, -809.3934203748788, ["cpx", "g", "opx", "ol"], [0.12366965263363426, 0.126211663721959, 0.13405802184212934, 0.6160606618022773]), outP{Float64}(30.0, 1400.0, 0, -835.0206475259749, ["cpx", "opx", "g", "ol"], [0.1727768966756921, 0.11084299123221208, 0.10031945028981894, 0.6160606618022767]), outP{Float64}(30.0, 1600.0, 0, -862.4847644122326, ["liq", "opx", "ol"], [0.21001279048221477, 0.24319492983069038, 0.546793222111941]), outP{Float64}(30.0, 1800.0, 0, -892.4957920942894, ["liq", "ol"], [0.6127489985757362, 0.38722526310731425]), outP{Float64}(30.0, 2000.0, 0, -925.6796244882266, ["liq"], [0.9999900671381493]), outP{Float64}(40.0, 800.0, 0, -748.5532733319261, ["g", "cpx", "opx", "ol"], [0.14739693263236445, 0.11096742588885528, 0.12557497806772447, 0.6160606618022769]), outP{Float64}(40.0, 1000.0, 0, -770.2508590133489, ["ol", "cpx", "opx", "g"], [0.6160606618022769, 0.11121580852269891, 0.12760561315561728, 0.14511791651940673]), outP{Float64}(40.0, 1200.0, 0, -793.9642970278298, ["opx", "g", "ol", "cpx"], [0.12644438133435348, 0.14088882800800623, 0.616060661802277, 0.11660612885536334]), outP{Float64}(40.0, 1400.0, 0, -819.4594227524026, ["opx", "g", "ol", "cpx"], [0.11141553948794847, 0.13199773100823958, 0.616060661802277, 0.14052606770153478]), outP{Float64}(40.0, 1600.0, 0, -846.5666746014005, ["opx", "cpx", "ol", "g"], [0.010134230286265417, 0.26979032137893766, 0.6160606618022769, 0.10401478653252011]), outP{Float64}(40.0, 1800.0, 0, -876.0425087378393, ["ol", "opx", "liq"], [0.39031262219050544, 0.13793131060477412, 0.4717232861072306]), outP{Float64}(40.0, 2000.0, 0, -908.8263682835493, ["liq"], [0.9999886578287357]), outP{Float64}(50.0, 800.0, 0, -733.4676382458092, ["ol", "cpx", "opx", "g"], [0.6160606618022768, 0.10957162568711479, 0.1254740005345372, 0.1488937119760829]), outP{Float64}(50.0, 1000.0, 0, -755.063428780673, ["opx", "ol", "g", "cpx"], [0.12767850038071737, 0.6160606618022773, 0.14834369035764555, 0.10791714745935985]), outP{Float64}(50.0, 1200.0, 0, -778.6687362464927, ["opx", "cpx", "g", "ol"], [0.12676123837083356, 0.11045655602601663, 0.1467215438008727, 0.616060661802277]), outP{Float64}(50.0, 1400.0, 0, -804.0468747476191, ["opx", "g", "ol", "cpx"], [0.11667250455079639, 0.14266440972411493, 0.616060661802277, 0.12460242392281141]), outP{Float64}(50.0, 1600.0, 0, -831.0179644511526, ["opx", "cpx", "g", "ol"], [0.05622968697275506, 0.19676053825872364, 0.13094911296624398, 0.6160606618022771]), outP{Float64}(50.0, 1800.0, 0, -859.8989781579687, ["ol", "liq", "opx"], [0.40889276245276723, 0.3873550753518144, 0.20372996376395922]), outP{Float64}(50.0, 2000.0, 0, -892.2515301871125, ["liq"], [0.9999877403386883])]

test1list       = Vector{outP}(undef, 42)
test1list       = outP[outP{Float64}(0.0, 800.0, 1, -943.3063927482982, ["ol", "pl4T", "cpx", "spn", "opx"], [0.056307003240191396, 0.4684458442389433, 0.3103597316381442, 0.0123452940574957, 0.15254212682522508]), outP{Float64}(0.0, 1000.0, 1, -969.1906195977808, ["ol", "pl4T", "cpx", "spn", "opx"], [0.054257918246936224, 0.46078528298226046, 0.34628892976102654, 0.01026946145099861, 0.12839840755877802]), outP{Float64}(0.0, 1200.0, 1, -998.1073123693945, ["liq", "ol", "pl4T"], [0.8427280466045326, 0.058798597367722784, 0.09849573305079129]), outP{Float64}(0.0, 1400.0, 1, -1031.9865613209831, ["liq"], [0.9999369931286983]), outP{Float64}(0.0, 1600.0, 1, -1068.192380824948, ["liq"], [0.9999652651386562]), outP{Float64}(0.0, 1800.0, 1, -1106.3970408669986, ["liq"], [0.9999652695461332]), outP{Float64}(0.0, 2000.0, 1, -1146.3926528475686, ["liq"], [0.9999577042177217]), outP{Float64}(10.0, 800.0, 1, -923.4514344178708, ["g", "cpx", "pl4T", "opx", "ru"], [0.15956382939360217, 0.3806172155361415, 0.32672483843750494, 0.12981168226769169, 0.0032824343650595346]), outP{Float64}(10.0, 1000.0, 1, -949.1596595854129, ["opx", "g", "cpx", "pl4T"], [0.16671351013689048, 0.02796212566940595, 0.4643756435076068, 0.34094872068609655]), outP{Float64}(10.0, 1200.0, 1, -977.3139093894723, ["cpx", "pl4T", "opx"], [0.6041719361217365, 0.3409487206860963, 0.054879343192167036]), outP{Float64}(10.0, 1400.0, 1, -1009.99773404958, ["liq"], [0.9999418408143267]), outP{Float64}(10.0, 1600.0, 1, -1045.9440455705483, ["liq"], [0.9999650109725381]), outP{Float64}(10.0, 1800.0, 1, -1083.8882408130862, ["liq"], [0.999960371621147]), outP{Float64}(10.0, 2000.0, 1, -1123.6211320023, ["liq"], [0.99996911316929]), outP{Float64}(20.0, 800.0, 1, -905.3450770225846, ["cpx", "g", "q", "ky", "ru"], [0.4735946660302741, 0.40401560679314247, 0.04994533392031061, 0.0699348050753592, 0.002509588180913409]), outP{Float64}(20.0, 1000.0, 1, -930.4514878864085, ["cpx", "g", "pl4T", "q", "ru"], [0.3971722305441107, 0.514370973514709, 0.012709218517705458, 0.07445084344392047, 0.0012967339795543432]), outP{Float64}(20.0, 1200.0, 1, -958.044788234852, ["g", "cpx", "pl4T", "q"], [0.35650002611622994, 0.4762831633885368, 0.11509723743797408, 0.05211957305725895]), outP{Float64}(20.0, 1400.0, 1, -988.8278413916834, ["liq", "cpx"], [0.7746364768039748, 0.2253878979558529]), outP{Float64}(20.0, 1600.0, 1, -1024.439613519149, ["liq"], [0.9999708939567556]), outP{Float64}(20.0, 1800.0, 1, -1062.1482427981314, ["liq"], [0.9999600688617434]), outP{Float64}(20.0, 2000.0, 1, -1101.6447802426374, ["liq"], [0.999971159101672]), outP{Float64}(30.0, 800.0, 1, -887.8693098936317, ["cpx", "g", "coe", "ky", "ru"], [0.38200415902357476, 0.5241578952143852, 0.0673543266087381, 0.024251954812089867, 0.0022316643412119996]), outP{Float64}(30.0, 1000.0, 1, -912.8505553436257, ["g", "cpx", "q", "ru"], [0.582980704068198, 0.33833882192731807, 0.07728151302665906, 0.0013989609778247627]), outP{Float64}(30.0, 1200.0, 1, -940.1448990398449, ["g", "cpx", "q"], [0.5464376628809569, 0.3748818631145593, 0.07868047400448375]), outP{Float64}(30.0, 1400.0, 1, -969.5877048388502, ["g", "liq", "cpx"], [0.40468311436618903, 0.22272434970860602, 0.3725900663018206]), outP{Float64}(30.0, 1600.0, 1, -1003.4856223541677, ["liq"], [0.999976658968688]), outP{Float64}(30.0, 1800.0, 1, -1040.973488906311, ["liq"], [0.9999683497926003]), outP{Float64}(30.0, 2000.0, 1, -1080.249147496561, ["liq"], [0.9999703752244993]), outP{Float64}(40.0, 800.0, 1, -870.6612865872631, ["g", "cpx", "coe", "ru"], [0.5899332096745897, 0.33138631632092624, 0.07717594600734427, 0.0015045279971395702]), outP{Float64}(40.0, 1000.0, 1, -895.5413701572314, ["g", "cpx", "coe", "ru"], [0.595765192151152, 0.325554333844364, 0.07813948661584325, 0.000540987388640556]), outP{Float64}(40.0, 1200.0, 1, -922.6855215570409, ["g", "cpx", "coe"], [0.5923350403628342, 0.32898448563268173, 0.07868047400448376]), outP{Float64}(40.0, 1400.0, 1, -951.8302281104818, ["cpx", "g", "coe"], [0.3645632245677071, 0.556756301427809, 0.07868047400448383]), outP{Float64}(40.0, 1600.0, 1, -983.3365803597649, ["cpx", "g", "liq"], [0.1939496297596014, 0.3214224034829167, 0.4846033858892608]), outP{Float64}(40.0, 1800.0, 1, -1020.2479864714834, ["liq"], [0.9999657672193685]), outP{Float64}(40.0, 2000.0, 1, -1059.3133154436234, ["liq"], [0.999963362778731]), outP{Float64}(50.0, 800.0, 1, -853.5780924421651, ["g", "cpx", "coe", "ru"], [0.5967949877596213, 0.3245245382358948, 0.07820521683727472, 0.00047525716720898926]), outP{Float64}(50.0, 1000.0, 1, -878.3640005086056, ["cpx", "g", "coe"], [0.32151734327869025, 0.5998021827168257, 0.07868047400448373]), outP{Float64}(50.0, 1200.0, 1, -905.4050507940307, ["cpx", "g", "coe"], [0.3221848308469694, 0.5991346951485466, 0.07868047400448376]), outP{Float64}(50.0, 1400.0, 1, -934.4268026304853, ["cpx", "g", "coe"], [0.3293460444602249, 0.5919734815352912, 0.07868047400448379]), outP{Float64}(50.0, 1600.0, 1, -965.2275762944403, ["cpx", "g", "coe"], [0.363255608425217, 0.5580639175702992, 0.07868047400448379]), outP{Float64}(50.0, 1800.0, 1, -999.8956828719806, ["liq"], [0.9999665411574561]), outP{Float64}(50.0, 2000.0, 1, -1038.7584782945137, ["liq"], [0.9999619565011438])]

test5list       = Vector{outP}(undef, 96)
test5list       = outP[outP{Float64}(0.0, 800.0, 5, -816.8408103758256, ["liq", "spn", "pl4T", "sph"], [0.8862912362779436, 0.014358927139203049, 0.09503222740459039, 0.0002734944916813597]), outP{Float64}(0.0, 1000.0, 5, -848.5007387121657, ["pl4T", "spn", "liq"], [0.004986299607829192, 0.004684011958261555, 0.9900071198012075]), outP{Float64}(0.0, 1200.0, 5, -883.1864513009938, ["liq"], [0.9999022073710113]), outP{Float64}(0.0, 1400.0, 5, -920.3900870879809, ["liq"], [1.0000401678861603]), outP{Float64}(0.0, 1600.0, 5, -959.6881327405104, ["liq"], [1.000049541931036]), outP{Float64}(0.0, 1800.0, 5, -1000.7960184839209, ["liq"], [0.9999945579314595]), outP{Float64}(0.0, 2000.0, 5, -1043.5201881257549, ["liq"], [1.0000062479854939]), outP{Float64}(10.0, 800.0, 5, -794.1719012103729, ["liq", "hb", "pl4T", "spn"], [0.8493560630485644, 0.06866428302736949, 0.07429331625775014, 0.007726425105909219]), outP{Float64}(10.0, 1000.0, 5, -825.1092138464046, ["liq"], [0.9999132784456582]), outP{Float64}(10.0, 1200.0, 5, -858.9339178710412, ["liq"], [0.9999978692854135]), outP{Float64}(10.0, 1400.0, 5, -895.0612865090908, ["liq"], [0.9999444858429375]), outP{Float64}(10.0, 1600.0, 5, -933.1032527466884, ["liq"], [1.000069213289366]), outP{Float64}(10.0, 1800.0, 5, -972.8043647705322, ["liq"], [1.0000211448162872]), outP{Float64}(10.0, 2000.0, 5, -1013.9733835801312, ["liq"], [1.0000269560744672]), outP{Float64}(20.0, 800.0, 5, -773.4423390882148, ["liq", "cpx", "g", "ep", "q", "ky"], [0.7717337023006312, 0.05019802754361492, 0.06651968206010657, 0.07250801582066095, 0.013230646877610814, 0.02580015301877901]), outP{Float64}(20.0, 1000.0, 5, -803.498035603907, ["liq"], [0.9999098054315866]), outP{Float64}(20.0, 1200.0, 5, -836.579964917351, ["liq"], [0.9999943557214976]), outP{Float64}(20.0, 1400.0, 5, -871.7992681864683, ["liq"], [1.0000546208430006]), outP{Float64}(20.0, 1600.0, 5, -908.8188247089591, ["liq"], [1.0000538016148903]), outP{Float64}(20.0, 1800.0, 5, -947.4279425120193, ["liq", "fl"], [0.9290133084854321, 0.07099012608673615]), outP{Float64}(20.0, 2000.0, 5, -987.6123267430705, ["liq", "fl"], [0.7979396634404049, 0.20206473889473647]), outP{Float64}(30.0, 800.0, 5, -754.303751505519, ["ep", "cpx", "g", "liq", "fl", "coe", "ky"], [0.02240180420665831, 0.23015312808540078, 0.019381542639129185, 0.4404280593203566, 0.1085220515932453, 0.09538788495747282, 0.083713403461023]), outP{Float64}(30.0, 1000.0, 5, -783.0412308339419, ["g", "liq", "cpx", "ky"], [0.020883926997142848, 0.9398668280246678, 0.02028657521096909, 0.01895361947346346]), outP{Float64}(30.0, 1200.0, 5, -815.3866836631354, ["liq"], [1.0000313107543697]), outP{Float64}(30.0, 1400.0, 5, -849.7848904225832, ["liq", "fl"], [0.9715512134992117, 0.02844723275430431]), outP{Float64}(30.0, 1600.0, 5, -886.0921189170638, ["liq", "fl"], [0.8173516328439473, 0.18265466008942627]), outP{Float64}(30.0, 1800.0, 5, -924.2309718221084, ["liq", "fl"], [0.6644581537688511, 0.3355458608818084]), outP{Float64}(30.0, 2000.0, 5, -964.0734651934788, ["liq", "fl"], [0.5232426277880491, 0.47675656343992123]), outP{Float64}(40.0, 800.0, 5, -736.599442189669, ["cpx", "g", "liq", "fl", "coe", "ky"], [0.2623351529185757, 0.06342671351515282, 0.1417686279782993, 0.29662939052313353, 0.14654374244463508, 0.08929519973684712]), outP{Float64}(40.0, 1000.0, 5, -763.9422896785777, ["g", "liq", "cpx", "fl", "coe", "ky"], [0.0035449299115794686, 0.4048797210253103, 0.20855606354475972, 0.19903189937141388, 0.11099184095314296, 0.07294976933125574]), outP{Float64}(40.0, 1200.0, 5, -795.0361989994902, ["liq", "fl"], [0.9649803886212627, 0.03496183468615838]), outP{Float64}(40.0, 1400.0, 5, -828.9164657806044, ["fl", "liq"], [0.20867202070367544, 0.7913282020049616]), outP{Float64}(40.0, 1600.0, 5, -864.9033972682595, ["liq", "fl"], [0.613481233551065, 0.38652091993862553]), outP{Float64}(40.0, 1800.0, 5, -902.8152892852502, ["liq", "fl"], [0.46111608368465395, 0.538880476523731]), outP{Float64}(40.0, 2000.0, 5, -942.4518074133127, ["liq", "fl"], [0.35224040122151024, 0.6477598197708965]), outP{Float64}(50.0, 800.0, 5, -719.5098791252326, ["cpx", "g", "liq", "fl", "coe", "ky"], [0.23537286660011253, 0.06667149282100356, 0.24004854427745212, 0.22443385729284537, 0.1555803094598784, 0.07787789740118886]), outP{Float64}(50.0, 1000.0, 5, -746.3435522442187, ["cpx", "g", "liq", "fl", "coe", "ky"], [0.24710051717231796, 0.06628434597452318, 0.1121838474941778, 0.3471756827631059, 0.14794239193704872, 0.07930897618385349]), outP{Float64}(50.0, 1200.0, 5, -775.9812400431259, ["cpx", "liq", "fl", "coe", "ky"], [0.1623028650591822, 0.38277436299972345, 0.29490570394094856, 0.11623462016964657, 0.04375004849886859]), outP{Float64}(50.0, 1400.0, 5, -809.0639445913106, ["liq", "fl", "coe"], [0.6016824648793869, 0.3669081294107678, 0.031408649666256096]), outP{Float64}(50.0, 1600.0, 5, -844.8397028922531, ["liq", "fl"], [0.4548394821378891, 0.5451577917619518]), outP{Float64}(50.0, 1800.0, 5, -882.5756914546351, ["liq", "fl"], [0.33267225148247115, 0.6673279943890249]), outP{Float64}(50.0, 2000.0, 5, -922.0115667449019, ["liq", "fl"], [0.2591064729290612, 0.7408932404195568]), outP{Float64}(0.0, 800.0, 5, -816.8408103758256, ["liq", "spn", "pl4T", "sph"], [0.8862912362779436, 0.014358927139203049, 0.09503222740459039, 0.0002734944916813597]), outP{Float64}(0.0, 850.0, 5, -824.671712369329, ["liq", "spn", "pl4T"], [1.0003757317486524, 0.010868538069772858, 0.0974159164496503]), outP{Float64}(0.0, 900.0, 5, -832.1909492349512, ["liq", "spn", "pl4T"], [0.9274160695802545, 0.011237269810223901, 0.061256561612430066]), outP{Float64}(0.0, 950.0, 5, -840.2860990322195, ["liq", "spn", "pl4T"], [1.0163878927284442, 0.008031782483379489, 0.03013676233939497]), outP{Float64}(0.0, 1000.0, 5, -848.5007387121657, ["pl4T", "spn", "liq"], [0.004986299607829192, 0.004684011958261555, 0.9900071198012075]), outP{Float64}(0.0, 1050.0, 5, -856.9087376712222, ["liq", "spn"], [0.9997239033981656, 0.0002085373537145206]), outP{Float64}(0.0, 1100.0, 5, -865.4986681473456, ["liq"], [0.9998592470032165]), outP{Float64}(0.0, 1150.0, 5, -874.2594682411669, ["liq"], [0.9999793018967215]), outP{Float64}(0.0, 1200.0, 5, -883.1864513009938, ["liq"], [0.9999022073710113]), outP{Float64}(5.0, 800.0, 5, -805.3015681432512, ["liq", "hb", "pl4T", "spn", "fl"], [0.6499026138326842, 0.08294176358380062, 0.12881150964300248, 0.008669435143927234, 0.1296943784132022]), outP{Float64}(5.0, 850.0, 5, -812.7333113113299, ["liq", "pl4T", "hb", "spn", "fl"], [0.7572416956171051, 0.11175961544874358, 0.027645374850337988, 0.010815109218688149, 0.09251321361345562]), outP{Float64}(5.0, 900.0, 5, -820.4347503199606, ["liq", "pl4T", "spn", "fl"], [0.8440510713789308, 0.08496321816120869, 0.010433688503885577, 0.06054331421016162]), outP{Float64}(5.0, 950.0, 5, -828.3666775893157, ["liq", "pl4T", "spn", "fl"], [0.9136169566117738, 0.046326255221672726, 0.007387110729204516, 0.03266324992719952]), outP{Float64}(5.0, 1000.0, 5, -836.5158442499047, ["liq", "spn", "fl"], [0.9964606635139137, 0.002408519699498244, 0.001070594684923099]), outP{Float64}(5.0, 1050.0, 5, -844.8176950543566, ["liq", "fl"], [0.9960992411658615, 0.0038327485722609803]), outP{Float64}(5.0, 1100.0, 5, -853.2960850464276, ["liq", "fl"], [0.9931576575186539, 0.006800612014895942]), outP{Float64}(5.0, 1150.0, 5, -861.9412994290075, ["liq", "fl"], [0.9899551478989997, 0.01006018278990125]), outP{Float64}(5.0, 1200.0, 5, -870.7440438601182, ["liq", "fl"], [0.986698904693493, 0.013315166925602832]), outP{Float64}(10.0, 800.0, 5, -794.1719012103729, ["liq", "hb", "pl4T", "spn"], [0.8493560630485644, 0.06866428302736949, 0.07429331625775014, 0.007726425105909219]), outP{Float64}(10.0, 850.0, 5, -801.5537333036225, ["liq", "pl4T", "spn", "hb"], [0.9000621149228123, 0.07202792120498794, 0.009601019878407892, 0.018330418933027427]), outP{Float64}(10.0, 900.0, 5, -809.211829922295, ["liq", "pl4T", "spn"], [0.961300995581437, 0.051060432822722504, 0.008127517567886568]), outP{Float64}(10.0, 950.0, 5, -817.0661549315, ["liq", "pl4T", "spn"], [0.9834253697737497, 0.012272549666740478, 0.004249555603266666]), outP{Float64}(10.0, 1000.0, 5, -825.1092138464046, ["liq"], [0.9999132784456582]), outP{Float64}(10.0, 1050.0, 5, -833.3207756858991, ["liq"], [0.9999367542559396]), outP{Float64}(10.0, 1100.0, 5, -841.701863011494, ["liq"], [0.9999106706860602]), outP{Float64}(10.0, 1150.0, 5, -850.2424258215583, ["liq"], [1.0000206180325049]), outP{Float64}(10.0, 1200.0, 5, -858.9339178710412, ["liq"], [0.9999978692854135]), outP{Float64}(15.0, 800.0, 5, -783.5916320384399, ["liq", "ep", "g", "spn", "hb", "ky"], [0.834053603741585, 0.08330390545622761, 0.051161796543503314, 0.00023332543513059908, 0.02976732093814969, 0.001378884463761559]), outP{Float64}(15.0, 850.0, 5, -790.8643638042907, ["liq", "ep", "g", "spn"], [0.886141643683738, 0.06834501279793335, 0.04418110505900476, 0.001315256362943371]), outP{Float64}(15.0, 900.0, 5, -798.3851112048403, ["liq", "spn", "g", "ep"], [0.9416564593260035, 0.0017162402137684558, 0.008330461466592899, 0.04834319574594434]), outP{Float64}(15.0, 950.0, 5, -806.1759022483775, ["liq", "ep", "spn", "opx"], [0.9968168669287523, 0.0012933375788154143, 0.0003954120556132974, 0.0014356961080946046]), outP{Float64}(15.0, 1000.0, 5, -814.13351178708, ["liq"], [0.9999775389290733]), outP{Float64}(15.0, 1050.0, 5, -822.2615015743439, ["liq"], [0.999954493041535]), outP{Float64}(15.0, 1100.0, 5, -830.5513862832328, ["liq"], [0.9999327769394369]), outP{Float64}(15.0, 1150.0, 5, -838.9943486605599, ["liq"], [0.999972115690674]), outP{Float64}(15.0, 1200.0, 5, -847.5818791371571, ["liq"], [0.999920604514048]), outP{Float64}(20.0, 800.0, 5, -773.4423390882148, ["liq", "cpx", "g", "ep", "q", "ky"], [0.7717337023006312, 0.05019802754361492, 0.06651968206010657, 0.07250801582066095, 0.013230646877610814, 0.02580015301877901]), outP{Float64}(20.0, 850.0, 5, -780.5916527986406, ["g", "liq", "ep", "ky"], [0.0737172325251027, 0.8472720406347934, 0.07203499461604138, 0.006982976436930247]), outP{Float64}(20.0, 900.0, 5, -787.981829529926, ["liq", "g", "ep", "ky"], [0.8946062176370505, 0.048444056085985086, 0.055250621968014293, 0.001660847610632045]), outP{Float64}(20.0, 950.0, 5, -795.6145381100428, ["liq", "g", "ep"], [0.9545525888095877, 0.013815318501934822, 0.03163659990256949]), outP{Float64}(20.0, 1000.0, 5, -803.498035603907, ["liq"], [0.9999098054315866]), outP{Float64}(20.0, 1050.0, 5, -811.5435774111891, ["liq"], [0.9998968589971392]), outP{Float64}(20.0, 1100.0, 5, -819.7446420142709, ["liq"], [0.9999352107850018]), outP{Float64}(20.0, 1150.0, 5, -828.0925629313261, ["liq"], [1.0000473436675625]), outP{Float64}(20.0, 1200.0, 5, -836.579964917351, ["liq"], [0.9999943557214976]), outP{Float64}(25.0, 800.0, 5, -763.6955977549014, ["g", "cpx", "ep", "liq", "q", "ky"], [0.03611320448472103, 0.14259532083062246, 0.04222148788679202, 0.6703418365434568, 0.049429568892482935, 0.05927565166917306]), outP{Float64}(25.0, 850.0, 5, -770.657152902155, ["g", "liq", "ep", "cpx", "q", "ky"], [0.0447719567750754, 0.7528230322907074, 0.04137562447208241, 0.09297163313784372, 0.023509801149287423, 0.044573672986422126]), outP{Float64}(25.0, 900.0, 5, -777.8874888410147, ["g", "ep", "liq", "cpx", "ky"], [0.053440232274798015, 0.03697085867014797, 0.8457697938957073, 0.03784103903766265, 0.0260049286311109]), outP{Float64}(25.0, 950.0, 5, -785.3865564738696, ["liq", "g", "ep", "ky"], [0.9188222928895089, 0.046507360729575176, 0.025679636506612268, 0.008980800902641311]), outP{Float64}(25.0, 1000.0, 5, -793.1396821868473, ["ep", "liq", "g"], [0.00533103826281062, 0.9868575224412784, 0.007754442097056966]), outP{Float64}(25.0, 1050.0, 5, -801.1042496873441, ["liq"], [0.999851086583135]), outP{Float64}(25.0, 1100.0, 5, -809.2174903440497, ["liq"], [0.9999089163621854]), outP{Float64}(25.0, 1150.0, 5, -817.4728406717431, ["liq"], [0.999929821330627]), outP{Float64}(25.0, 1200.0, 5, -825.8632965653837, ["liq"], [0.9999697639727972])]

gv.verbose  = -1    # switch off any verbose
@testset verbose = true "Total tests" begin
    @testset "KLB-1 peridorite tests" begin
        for i=1:size(test0list,1)
            
            bulk_rock   = get_bulk_rock(gv, test0list[i].test)
            out         = point_wise_minimization(test0list[i].P,test0list[i].T, bulk_rock, gv, DB)

            @test out.G_system  ≈ test0list[i].G
            @test out.ph        == test0list[i].ph
            @test all(abs.(out.ph_frac - test0list[i].ph_frac)  .< 1e-6)

            @show out
        end
    end

    @testset "RE-46 icelandic basalt tests" begin
        for i=1:size(test1list,1)
            
            bulk_rock   = get_bulk_rock(gv, test1list[i].test)
            out         = point_wise_minimization(test1list[i].P,test1list[i].T, bulk_rock, gv, DB)

            @test out.G_system  ≈ test1list[i].G
            @test out.ph        == test1list[i].ph
            @test all(abs.(out.ph_frac - test1list[i].ph_frac)  .< 1e-6)

            @show out
        end
    end

    @testset "WET MORB tests" begin
        for i=1:size(test5list,1)
            
            bulk_rock   = get_bulk_rock(gv, test5list[i].test)
            out         = point_wise_minimization(test5list[i].P,test5list[i].T, bulk_rock, gv, DB)

            @test out.G_system  ≈ test5list[i].G
            @test out.ph        == test5list[i].ph
            @test all(abs.(out.ph_frac - test5list[i].ph_frac)  .< 1e-6)

            @show out
        end
    end
end
