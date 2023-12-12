%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test script for the ADG extrapolation model code. The ADG ext code is run
%for one specified input of ag and adg at specified hyperspectral
%wavelengths (lambda). The resulting output from the test script is saved
%to ADG_ext_test_run_yyyymmdd.xls for comparison with the provided output
%file ADG_ext_test_run.xls.
%
%Reference:
%
%Kehrli, M. D., Stramski, D., Reynolds, R. A., & Joshi, I. D. (2023).
%Estimation of chromophoric dissolved organic matter and non-algal
%particulate absorption coefficients of seawater in the ultraviolet by
%extrapolation from the visible spectral region. Optics Express, 31(11),
%17450. https://doi.org/10.1364/OE.486354
%
%Created: July 24, 2023
%Completed: July 31, 2023
%Updates:
%
%M. D. Kehrli, D. Stramski, R. A. Reynolds, and I. D. Joshi
%Ocean Optics Research Laboratory, Scripps Institution of Oceanography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear command window and workspace; close figures
clc; clearvars; close all;

%define input parameters

%input ag [m^-1]
ag = [0.152207861000000	0.149720621000000 0.147211790000000...
    0.144668415000000 0.142113524000000	0.139490983000000...
    0.137629504000000 0.135418621000000	0.133243255000000...
    0.131102833000000 0.128996795000000	0.126924588000000...
    0.124885669000000 0.122879504000000	0.120905566000000...
    0.118963337000000 0.117052308000000	0.115171978000000...
    0.113321853000000 0.111501449000000	0.109710288000000...
    0.107947900000000 0.106213823000000	0.104507603000000...
    0.102828791000000 0.101176948000000	0.0995516390000000...
    0.0979524400000000 0.0963789310000000 0.0948306980000000...
    0.0933073360000000 0.0918084460000000 0.0903336330000000...
    0.0888825120000000 0.0874547020000000 0.0860498280000000...
    0.0846675230000000 0.0833074220000000 0.0819691700000000...
    0.0806524160000000 0.0793568150000000 0.0780820250000000...
    0.0768277140000000 0.0755935530000000 0.0743792170000000...
    0.0731843880000000 0.0720087530000000 0.0708520030000000...
    0.0697138350000000 0.0685939510000000 0.0674920570000000...
    0.0664078630000000 0.0653410860000000 0.0642914460000000...
    0.0632586670000000 0.0622424790000000 0.0612426150000000...
    0.0602588120000000 0.0592908140000000 0.0583383650000000...
    0.0574012170000000 0.0564791230000000 0.0555718410000000...
    0.0546791340000000 0.0538007680000000 0.0529365110000000...
    0.0520861380000000 0.0512494260000000 0.0504261540000000...
    0.0496161080000000 0.0488190740000000 0.0480348430000000...
    0.0472632110000000 0.0465039740000000 0.0457569330000000...
    0.0450218930000000 0.0442986610000000 0.0435870460000000...
    0.0428868630000000 0.0421979280000000 0.0415200600000000...
    0.0408530810000000 0.0401968160000000 0.0395510940000000...
    0.0389157450000000 0.0382906010000000 0.0376755010000000...
    0.0370702810000000 0.0364747830000000 0.0358888520000000...
    0.0353123320000000 0.0347450750000000 0.0341869290000000...
    0.0336377500000000 0.0330973920000000 0.0325657150000000...
    0.0320425790000000 0.0315278470000000 0.0310213830000000...
    0.0305230550000000 0.0300327320000000 0.0295502860000000...
    0.0290755890000000 0.0286085190000000 0.0281489510000000...
    0.0276967660000000 0.0272518450000000 0.0268140700000000...
    0.0263833290000000 0.0259595070000000 0.0255424930000000...
    0.0251321780000000 0.0247284540000000 0.0243312160000000...
    0.0239403590000000 0.0235557800000000 0.0231773800000000...
    0.0228050580000000 0.0224387170000000 0.0220782610000000...
    0.0217235960000000 0.0213746280000000 0.0210312650000000...
    0.0206934190000000 0.0203609990000000 0.0200339200000000...
    0.0197120950000000 0.0193954390000000 0.0190838710000000...
    0.0187773070000000 0.0184756680000000 0.0181788750000000...
    0.0178868490000000 0.0175995140000000 0.0173167950000000...
    0.0170386180000000 0.0167649090000000 0.0164955980000000...
    0.0162306120000000 0.0159698830000000 0.0157133430000000...
    0.0154609230000000 0.0152125590000000 0.0149681840000000...
    0.0147277350000000 0.0144911480000000 0.0142583620000000...
    0.0140293160000000 0.0138039490000000 0.0135822020000000...
    0.0133640170000000 0.0131493370000000 0.0129381060000000...
    0.0127302680000000 0.0125257690000000 0.0123245550000000...
    0.0121265730000000 0.0119317710000000 0.0117400990000000...
    0.0115515060000000 0.0113659430000000 0.0111833600000000...
    0.0110037100000000 0.0108269460000000 0.0106530220000000...
    0.0104818920000000 0.0103135110000000 0.0101478340000000...
    0.00998481900000000 0.00982442300000000	0.00966660300000000...
    0.00951131900000000	0.00935852900000000	0.00920819300000000...
    0.00906027300000000	0.00891472800000000	0.00877152200000000...
    0.00863061600000000	0.00849197400000000	0.00835555800000000...
    0.00822133500000000	0.00808926700000000	0.00795932100000000...
    0.00783146200000000	0.00770565700000000	0.00758187400000000...
    0.00746007800000000	0.00734024000000000	0.00722232600000000...
    0.00710630600000000	0.00699215000000000	0.00687982800000000...
    0.00676931100000000	0.00666056800000000	0.00655357300000000...
    0.00644829600000000	0.00634471100000000	0.00624278900000000...
    0.00614250500000000	0.00604383200000000	0.00594674300000000...
    0.00585121500000000	0.00575722100000000	0.00566473700000000...
    0.00557373800000000	0.00548420200000000	0.00539610300000000...
    0.00530942000000000	0.00522413000000000	0.00514020900000000...
    0.00505763700000000	0.00497639100000000	0.00489645000000000...
    0.00481779300000000	0.00474040000000000	0.00466425000000000...
    0.00458932400000000	0.00451560100000000	0.00444306200000000...
    0.00437168900000000	0.00430146200000000	0.00423236300000000...
    0.00416437400000000	0.00409747800000000	0.00403165600000000...
    0.00396689100000000	0.00390316700000000	0.00384046600000000...
    0.00377877300000000	0.00371807100000000	0.00365834400000000...
    0.00359957600000000	0.00354175200000000	0.00348485800000000...
    0.00342887700000000	0.00337379500000000	0.00331959900000000...
    0.00326627200000000	0.00321380300000000	0.00316217600000000...
    0.00311137900000000	0.00306139800000000	0.00301222000000000...
    0.00296383100000000	0.00291622000000000	0.00286937400000000...
    0.00282328000000000	0.00277792700000000	0.00273330200000000...
    0.00268939500000000	0.00264619200000000	0.00260368400000000...
    0.00256185800000000	0.00252070400000000	0.00248021200000000...
    0.00244037000000000	0.00240116700000000	0.00236259500000000...
    0.00232464200000000	0.00228729900000000	0.00225055600000000...
    0.00221440300000000	0.00217883100000000	0.00214383000000000...
    0.00210939100000000	0.00207550600000000	0.00204216500000000...
    0.00200936000000000	0.00197708100000000	0.00194532200000000...
    0.00191407200000000	0.00188332400000000	0.00185307000000000...
    0.00182330300000000	0.00179401300000000	0.00176519400000000...
    0.00173683800000000	0.00170893700000000	0.00168148500000000...
    0.00165447400000000	0.00162789600000000	0.00160174600000000...
    0.00157601500000000	0.00155069800000000	0.00152578800000000...
    0.00150127700000000	0.00147716100000000	0.00145343200000000...
    0.00143008400000000	0.00140711100000000	0.00138450700000000...
    0.00136226600000000	0.00134038300000000	0.00131885100000000...
    0.00129766500000000	0.00127681900000000	0.00125630800000000...
    0.00123612700000000	0.00121627000000000	0.00119673100000000...
    0.00117750700000000]';

%input adg [m^-1]
ad = [0.0397546087688120 0.0394517909845550	0.0391550533752260...
    0.0388627910574420 0.0385730362469460 0.0382840737833890...
    0.0379949408387830 0.0377060315822230 0.0374179173693950...
    0.0371311844247890 0.0368466334772670 0.0365651195559840...
    0.0362867104451970 0.0360108664218730 0.0357366530107940...
    0.0354630121944680 0.0351891716696860 0.0349141808533600...
    0.0346370412761000 0.0343572327484040 0.0340752306347010...
    0.0337916027192490 0.0335072127775580 0.0332234305618150...
    0.0329414723256630 0.0326621447017560 0.0323860586959260...
    0.0321127097163340 0.0318409331128350 0.0315697353723100...
    0.0312976248037970 0.0310236077484040 0.0307475089145840...
    0.0304693866842640 0.0301900320195410 0.0299109661303280...
    0.0296334058533600 0.0293589437542350 0.0290889522819310...
    0.0288244754597740 0.0285666555618150 0.0283163502411150...
    0.0280738632877620 0.0278396935355760 0.0276136337688120...
    0.0273952342061300 0.0271841319466544 0.0269794568008818...
    0.0267801947746428 0.0265857466696865 0.0263951120486952...
    0.0262075675880539 0.0260228458679373 0.0258402981274124...
    0.0256598549058381 0.0254818001682288 0.0253060930982579...
    0.0251331490749343 0.0249636893810568 0.0247977477629810...
    0.0246355046871792 0.0244769289583162 0.0243214538125437...
    0.0241687324568585 0.0240181201390743 0.0238686401828060...
    0.0237198087688119 0.0235714668591909 0.0234231419320772...
    0.0232750829670626 0.0231276013344095 0.0229810424422812...
    0.0228360723985495 0.0226933597892200 0.0225529922965086...
    0.0224152654743512 0.0222799278650218 0.0221462780836807...
    0.0220137573110859 0.0218814685355757 0.0217485036667710...
    0.0216143757513191 0.0214787122673541 0.0213411088416982...
    0.0212016352265378 0.0210601261886369 0.0209165904014649...
    0.0207708843519022 0.0206229279379081 0.0204727515530684...
    0.0203208696288702 0.0201678229087536 0.0200145049058381...
    0.0198619408387827 0.0197110221798906 0.0195624152556923...
    0.0194163811449052 0.0192728195559839 0.0191313803431559...
    0.0189915737105028 0.0188527743664795 0.0187148839145845...
    0.0185777979816398 0.0184416076755174 0.0183065124131267...
    0.0181726054160422 0.0180397328212900 0.0179076998037973...
    0.0177759738562754 0.0176442710137098 0.0175126186084620...
    0.0173809901828060 0.0172497512615232 0.0171195647454882...
    0.0169909917863046 0.0168643743664795 0.0167398859554008...
    0.0166172266259547 0.0164961533023396 0.0163761817279955...
    0.0162568021361588 0.0161376355909693 0.0160186395997157...
    0.0158997097163337 0.0157809570195407 0.0156623854451967...
    0.0155439907658964 0.0154256393081705 0.0153073301244970...
    0.0151890828212900 0.0150709957221646 0.0149532729087535...
    0.0148363382877623 0.0147205863198323 0.0146064053431559...
    0.0144940988562754 0.0143838398912608 0.0142757990749343...
    0.0141699405472375 0.0140661341332433 0.0139645296142929...
    0.0138654352994241 0.0137690843519022 0.0136756785938848...
    0.0135853348621063 0.0134981576026311 0.0134139437542346...
    0.0133321632148760 0.0132523080399489 0.0131738021361588...
    0.0130961118300364 0.0130188545414066 0.0129417350078789...
    0.0128648270632725 0.0127882569466544 0.0127117547600655...
    0.0126353688271209 0.0125591759699780 0.0124831096434474...
    0.0124072564364503 0.0123317165676457 0.0122565334043804...
    0.0121820808533600 0.0121083508970918 0.0120353652556923...
    0.0119632879233308 0.0118922082586078 0.0118219089145845...
    0.0117525421507361 0.0116842586959256 0.0116172827484037...
    0.0115514471798906 0.0114866039583162 0.0114222804160422...
    0.0113583650370334 0.0112944808533600 0.0112304640166253...
    0.0111663616113775 0.0111029134335349 0.0110406584772667...
    0.0109805161303279 0.0109231056347011 0.0108691350078789...
    0.0108187502411151 0.0107718325297448 0.0107277360282871...
    0.0106859756784329 0.0106459126317856 0.0106070965968002...
    0.0105688943373250 0.0105310632148760 0.0104938395997157...
    0.0104579035209985 0.0104236659845553 0.0103912795414066...
    0.0103606352994241 0.0103316030107944 0.0103036892352842...
    0.0102758252411151 0.0102464272819314 0.0102142801244970...
    0.0101788027192492 0.0101399016259547 0.0100978366113775...
    0.0100535011157506 0.0100080929524853 0.00996296328776230...
    0.00991928005161070	0.00987767837522590	0.00983825578047370...
    0.00980069083878270	0.00976416496414710	0.00972764710700430...
    0.00969043522653780	0.00965197618863690	0.00961210213615880...
    0.00957101671341830	0.00952922210700430	0.00948721314198970...
    0.00944571510991970	0.00940503785044450	0.00936535439563400...
    0.00932653624694600	0.00928837793790810	0.00925065818572150...
    0.00921354411866610	0.00917718260263110	0.00914199601370980...
    0.00910838690292260	0.00907680803994890	0.00904750709242700...
    0.00902044083878270	0.00899520665510920	0.00897127334607130...
    0.00894787203411800	0.00892424645102760	0.00889982669884100...
    0.00887432225277690	0.00884785760263110	0.00882061948309760...
    0.00879287035773310	0.00876516583878270	0.00873805840438040...
    0.00871165687376810	0.00868584003703340	0.00866007771924920...
    0.00863399805452610	0.00860713376881190	0.00857870301079440...
    0.00854799083878270	0.00851488967260200	0.00847946401662530...
    0.00844215410408880	0.00840335767551740	0.00836385745685850...
    0.00832480104286430	0.00828719105744160	0.00825168194665440...
    0.00821898245685850	0.00818963595540080	0.00816402422070690...
    0.00814203785044450	0.00812325191749990	0.00810754200496350...
    0.00809456919155240	0.00808346941021130	0.00807307903120250...
    0.00806211102828710	0.00804910592624630	0.00803262072216463...
    0.00801090760263110	0.00798286110117338	0.00794828814198970...
    0.00790769477464277	0.00786213231108591	0.00781368551808300...
    0.00776469156764568	0.00771757305452615	0.00767416357930749...
    0.00763562130525501	0.00760244652391390	0.00757449725277688...
    0.00755059193207717	0.00752931285044452	0.00750919666968650...
    0.00748897728193140	0.00746751838980312	0.00744399397289350...
    0.00741826314198970	0.00739040388542994	0.00736049776298096...
    0.00732883602828708	0.00729565031400137	0.00726120650933664...
    0.00722577582420545	0.00718903027026959	0.00715136540146492...
    0.00711330126152323	0.00707516343353489	0.00703752662595472...
    0.00700089448309757]';

%input wavelengths [nm]
lambda = (400:700)';

%input ADG extension correction LUT
load 'adg_cor_LUT.mat'

%extend ag, ad, and adg spectra using extension model
[lamout,agout,adout,adgout] = Extrapolate_ADG_VIStoUV(lambda,ag,ad,adg_cor_LUT);

%save inputs and outputs into an excel file
T1 = table(lambda,ag,ad);
T2 = table(lamout,agout,adout,adgout);
T1.Properties.VariableNames = {'Input Wavelength [nm]','Input ag [1/m]','Input ad [1/m]'};
T2.Properties.VariableNames = {'Output Wavelength [nm]','Output ag [1/m]','Output ad [1/m]','Output adg [1/m]'};
FormatOut = 'yyyymmdd';
outfile = [cd '\ADG_ext_test_run_' datestr(datetime,FormatOut)];
writetable(T1,outfile,'FileType','spreadsheet','Sheet','ADG_ext_Input')
writetable(T2,outfile,'FileType','spreadsheet','Sheet','ADG_ext_Output')