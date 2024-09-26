% Various commonly used lists

animalsNeuronexusUOL1 = {
  'M180313';
  'M180711_MD';
  'M180712_MD';
  'M180713_MD';
  'M180717_B_MD';
  'M180719_MD';
  'M180924_A_MD';
  'M180924_B_MD';
  'M181112_B_MD'};

animalsNeuronexusUOL2 = {
  'M181210_B_MD';
  'M190114_A_MD';
  'M190128_A_MD';
  'M190128_B_MD';
  'M190128_C_MD';
  'M190322_B_MD';
  'M190322_C_MD';
  'M190503_A_MD';
  'M190503_B_MD';
  'M190523_A_MD';
  'M190523_B_MD'};

animalsNeuropixelsUOL = {
  'M191018_MD';
  'M191106_MD';
  'M191107_MD';
  'M191119_A_MD';
  'M191119_B_MD';
  'M191119_C_MD';
  'M191128_A_MD';
  'M191128_B_MD';
  'M200316_MD';
  'M200317_MD';
  'M200318_MD';
  'M200319_MD';
  'M200323_MD';
  'M200324_MD';
  'M200325_MD';
  'M210802_MD'};

animalsAllensdk = {
  'M766640955';
  'M767871931';
  'M768515987';
  'M771160300';
  'M771990200';
  'M774875821';
  'M778240327';
  'M778998620';
  'M779839471';
  'M781842082';
  'M786091066';
  'M787025148';
  'M789848216';
  'M793224716';
  'M794812542';
  'M816200189';
  'M819186360';
  'M819701982';
  'M821695405';
  'M829720705';
  'M831882777';
  'M835479236';
  'M839068429';
  'M839557629';
  'M840012044';
  'M847657808'};

animalsNeuronexusUOL = [animalsNeuronexusUOL1; animalsNeuronexusUOL2];

animalsFull = [animalsNeuronexusUOL; animalsNeuropixelsUOL; animalsAllensdk];

animalsUOL = [animalsNeuronexusUOL; animalsNeuropixelsUOL];

animalsUOLReduced = [animalsNeuronexusUOL2; animalsNeuropixelsUOL];

animalsNeuronexusOI = animalsNeuronexusUOL;

animalsUOLOI = [animalsNeuronexusOI; animalsNeuropixelsUOL];

animalsNeuropixels = [animalsNeuropixelsUOL; animalsAllensdk];

animalsNeuropixelsOI = animalsNeuropixels;

animalsOI = [animalsNeuronexusOI; animalsNeuropixelsOI];

conditions = {
  'awake';
  'anaesthesia';
  'all'};

areas = {
  'lS1';
  'lVB1';
  'lPo';
  'lLP1';
  'lTh1';
  'lDG';
  'lCA1';
  'lHp1';
  'lRSC';
  'lVB2';
  'lLP2';
  'lLGN';
  'lTh2';
  'lCA3';
  'lVB';
  'lLP';
  'lTh';
  'lCA';
  'lHp';
  'lmPFC';
  'lV1';
  'rV1';
  'rVB';
  'rPo';
  'rLP';
  'rLGN'
  'rTh';
  'rDG';
  'rCA1';
  'rCA3';
  'rCA';
  'rHp';
  'VB';
  'LGN';
  'Th';
  'DG';
  'CA1';
  'CA3';
  'CA';
  'Hp';
  'V1';
  'Cx';
  'lVIS';
  'rVIS';
  'VIS';
  'S1';
  'RSC';
  'LP';
  'lV2';
  'rV2';
  'V2';
  'Po'};

areasOI = {
  'S1';
  'VB';
  'LGN';
  'LP';
  'Th';
  'DG';
  'RSC';
  'CA';
  'Hp';
  'V1';
  'V2';
  'VIS';
  'Po';
  'Cx'};
%areasOI = areas;

%iAreasOI = 1:numel(areas);
iAreasOI = find(ismember(areas,areasOI));

areasCritical = areasOI;

areasOIReduced = {
  'S1';
  'Th';
  'DG';
  'RSC';
  'CA';
  'Hp';
  'V1';
  'V2';
  'VIS';
  'Cx'};

areasCriticalReduced = areasOIReduced;

% areas2compare = {
%   'lVB1VslS1';
%   'lVB1VslDG';
%   'lVB1VslCA1';
%   'lVB1VslRSC';
%   'lVB1VslVB2';
%   'lVB1VslCA3';
%   'lVB1VslCA';
%   'lVB1VslHp';
%   'lVB2VslS1';
%   'lVB2VslDG';
%   'lVB2VslCA1';
%   'lVB2VslRSC';
%   'lVB2VslCA3';
%   'lVB2VslCA';
%   'lVB2VslHp';
%   'lVBVslS1';
%   'lVBVslDG';
%   'lVBVslCA1';
%   'lVBVslRSC';
%   'lVBVslCA3';
%   'lVBVslCA';
%   'lVBVslHp';
%   'lTh1VslS1';
%   'lTh1VslDG';
%   'lTh1VslCA1';
%   'lTh1VslRSC';
%   'lTh1VslTh2';
%   'lTh1VslCA3';
%   'lTh1VslCA';
%   'lTh1VslHp';
%   'lTh2VslS1';
%   'lTh2VslDG';
%   'lTh2VslCA1';
%   'lTh2VslRSC';
%   'lTh2VslCA3';
%   'lTh2VslCA';
%   'lTh2VslHp';
%   'lThVslS1';
%   'lThVslDG';
%   'lThVslCA1';
%   'lThVslRSC';
%   'lThVslTh2';
%   'lThVslCA3';
%   'lThVslCA';
%   'lThVslHp';
%   'lS1VslDG';
%   'lS1VslCA1';
%   'lS1VslRSC';
%   'lS1VslCA3';
%   'lS1VslCA';
%   'lS1VslHp';
%   'lDGVslCA1';
%   'lDGVslRSC';
%   'lDGVslCA3';
%   'lCA1VslRSC'
%   'lCA1VslCA3';
%   'lRSCVslCA3';
%   'lRSCVslCA';
%   'lRSCVslHp';
%   'lLGNVslV1';
%   'lLGNVslDG';
%   'lLGNVslCAl';
%   'lLGNVslCA3';
%   'lLGNVslCA';
%   'lLGNVslHp';
%   'lThVslV1';
%   'lV1VslDG';
%   'lV1VslCA1';
%   'lV1VslCA3';
%   'lV1VslCA';
%   'lV1VslHp';
%   'rVBVsrS1';
%   'rVBVsrDG';
%   'rVBVsrCA1';
%   'rVBVsrRSC';
%   'rVBVsrCA3';
%   'rVBVsrCA';
%   'rVBVsrHp';
%   'rThVsrS1';
%   'rThVsrDG';
%   'rThVsrCA1';
%   'rThVsrRSC';
%   'rThVsrTh2';
%   'rThVsrCA3';
%   'rThVsrCA';
%   'rThVsrHp';
%   'rS1VsrDG';
%   'rS1VsrCA1';
%   'rS1VsrRSC';
%   'rS1VsrCA3';
%   'rS1VsrCA';
%   'rS1VsrHp';
%   'rDGVsrCA1';
%   'rDGVsrRSC';
%   'rDGVsrCA3';
%   'rCA1VsrRSC'
%   'rCA1VsrCA3';
%   'rRSCVsrCA3';
%   'rRSCVsrCA';
%   'rRSCVsrHp';
%   'rLGNVsrV1';
%   'rLGNVsrDG';
%   'rLGNVsrCA1';
%   'rLGNVsrCA3';
%   'rLGNVsrCA';
%   'rLGNVsrHp';
%   'rThVsrV1';
%   'rV1VsrCA1';
%   'rV1VsrCA3';
%   'rV1VsrCA';
%   'rV1VsrHp';
%   'VBVsS1';
%   'VBVsDG';
%   'VBVsCA1';
%   'VBVsRSC';
%   'VBVsCA3';
%   'VBVsCA';
%   'VBVsHp';
%   'LGNVsV1';
%   'LGNVsDG';
%   'LGNVsCA1';
%   'LGNVsCA3';
%   'LGNVsCA';
%   'LGNVsHp';
%   'ThVsS1';
%   'ThVsDG';
%   'ThVsCA1';
%   'ThVsRSC';
%   'ThVsCA3';
%   'ThVsCA';
%   'ThVsHp';
%   'S1VsCA';
%   'S1VsHp';
%   'RSCVsCA';
%   'RSCVsHp'
%   'V1VsCA1';
%   'V1VsCA3';
%   'V1VsCA';
%   'V1VsHp';
%   'lDGVslCA';
%   'lCAVslDG';
%   'rDGVsrCA';
%   'rCAVsrDG';
%   'DGVsCA';
%   'CAVsDG';
%   'lDGVslV1';
%   'rDGVsrV1';
%   'DGVsV1';
%   'V1VsDG'};
areas2compare = {
  'VBVsS1';
  'VBVsRSC';
  'VBVsCA1';
  'VBVsCA3';
  'VBVsCA';
  'VBVsDG';
  'VBVsHp';
  'LGNVsV1';
  'LGNVsVIS';
  'LGNVsCA1';
  'LGNVsCA3';
  'LGNVsCA';
  'LGNVsDG';
  'LGNVsHp';
  'ThVsS1';
  'ThVsRSC';
  'ThVsV1';
  'ThVsVIS';
  'ThVsCA1';
  'ThVsCA3';
  'ThVsCA';
  'ThVsDG';
  'ThVsHp';
  'S1VsRSC';
  'S1VsCA1';
  'S1VsCA3';
  'S1VsCA';
  'S1VsDG';
  'S1VsHp';
  'RSCVsCA1';
  'RSCVsCA3';
  'RSCVsCA';
  'RSCVsDG';
  'RSCVsHp';
  'V1VsVIS';
  'V1VsCA1';
  'V1VsCA3';
  'V1VsCA';
  'V1VsDG';
  'V1VsHp';
  'VISVsCA1';
  'VISVsCA3';
  'VISVsCA';
  'VISVsDG';
  'VISVsHp';
  'CA1VsCA3';
  'CA1VsDG';
  'CA3VsDG';
  'CAVsDG';
  'S1VsVB';
  'RSCVsVB';
  'CA1VsVB';
  'CA3VsVB';
  'CAVsVB';
  'DGVsVB';
  'HpVsVB';
  'V1VsLGN';
  'VISVsLGN';
  'CA1VsLGN';
  'CA3VsLGN';
  'CAVsLGN';
  'DGVsLGN';
  'HpVsLGN';
  'S1VsTh';
  'RSCVsTh';
  'V1VsTh';
  'VISVsTh';
  'CA1VsTh';
  'CA3VsTh';
  'CAVsTh';
  'DGVsTh';
  'HpVsTh';
  'RSCVsS1';
  'CA1VsS1';
  'CA3VsS1';
  'CAVsS1';
  'DGVsS1';
  'HpVsS1';
  'CA1VsRSC';
  'CA3VsRSC';
  'CAVsRSC';
  'DGVsRSC';
  'HpVsRSC';
  'VISVsV1';
  'CA1VsV1';
  'CA3VsV1';
  'CAVsV1';
  'DGVsV1';
  'HpVsV1';
  'CA1VsVIS';
  'CA3VsVIS';
  'CAVsVIS';
  'DGVsVIS';
  'HpVsVIS';
  'CA3VsCA1';
  'DGVsCA1';
  'DGVsCA3';
  'DGVsCA';
  'LPVsLGN';
  'LPVsV1';
  'LPVsV2';
  'LPVsVIS';
  'LPVsCA1';
  'LPVsCA3';
  'LPVsCA';
  'LPVsDG';
  'LPVsHp';
  'LGNVsLP';
  'V1VsLP';
  'V2VsLP';
  'VISVsLP';
  'CA1VsLP';
  'CA3VsLP';
  'CAVsLP';
  'DGVsLP';
  'HpVsLP';
  'LGNVsV2';
  'ThVsV2';
  'V1VsV2';
  'V2VsV1';
  'V2VsVIS';
  'V2VsCA1';
  'V2VsCA3';
  'V2VsCA';
  'V2VsDG';
  'V2VsHp';
  'V2VsLGN';
  'V2VsTh';
  'VISVsV2';
  'CA1VsV2';
  'CA3VsV2';
  'CAVsV2';
  'DGVsV2';
  'HpVsV2';
  'PoVsVB';
  'PoVsS1';
  'PoVsRSC';
  'PoVsCA1';
  'PoVsCA3';
  'PoVsCA';
  'PoVsDG';
  'PoVsHp'
  'VBVsPo';
  'S1VsPo';
  'RSCVsPo';
  'CA1VsPo';
  'CA3VsPo';
  'CAVsPo';
  'DGVsPo';
  'HpVsPo';
  'CxVsTh';
  'CxVsHp';
  'ThVsCx';
  'HpVsCx'};

% areas2compareCritical = {
%   'VBVslS1';
%   'lS1VsVB';
%   'VBVslRSC';
%   'lRSCVsVB';
%   'VBVsCA';
%   'CAVsVB';
%   'lS1VslRSC';
%   'lRSCVslS1';
%   'lS1VsCA';
%   'CAVslS1';
%   'lRSCVsCA';
%   'CAVslRSC';
%   'VBVsDG';
%   'DGVsVB';
%   'VBVsCA1';
%   'CA1VsVB';
%   'VBVsCA3';
%   'CA3VsVB';
%   'LGNVsV1';
%   'V1VsLGN';
%   'LGNVsCA';
%   'CAVsLGN';
%   'V1VsCA';
%   'CAVsV1';
%   'ThVslS1';
%   'lS1VsTh';
%   'ThVslRSC';
%   'lRSCVsTh';
%   'ThVsCA';
%   'CAVsTh';
%   'ThVsDG';
%   'DGVsTh';
%   'ThVsCA1';
%   'CA1VsTh';
%   'ThVsCA3';
%   'CA3VsTh';
%   'ThVsV1';
%   'V1VsTh'};
areas2compareCritical = {
  'VBVsS1';
  'VBVsRSC';
  'VBVsCA';
  'VBVsDG';
  'VBVsHp';
  'LGNVsV1';
  'LGNVsVIS';
  'LGNVsCA';
  'LGNVsDG';
  'LGNVsHp';
  'ThVsS1';
  'ThVsRSC';
  'ThVsV1';
  'ThVsVIS';
  'ThVsCA';
  'ThVsDG';
  'ThVsHp';
  'S1VsRSC';
  'S1VsCA';
  'S1VsDG';
  'S1VsHp';
  'RSCVsCA';
  'RSCVsDG';
  'RSCVsHp';
  'V1VsVIS';
  'V1VsCA';
  'V1VsDG';
  'V1VsHp';
  'VISVsCA';
  'VISVsDG';
  'VISVsHp';
  'CAVsDG';
  'S1VsVB';
  'RSCVsVB';
  'CAVsVB';
  'DGVsVB';
  'HpVsVB';
  'V1VsLGN';
  'VISVsLGN';
  'CAVsLGN';
  'DGVsLGN';
  'HpVsLGN';
  'S1VsTh';
  'RSCVsTh';
  'V1VsTh';
  'VISVsTh';
  'CAVsTh';
  'DGVsTh';
  'HpVsTh';
  'RSCVsS1';
  'CAVsS1';
  'DGVsS1';
  'HpVsS1';
  'CAVsRSC';
  'DGVsRSC';
  'HpVsRSC';
  'VISVsV1';
  'CAVsV1';
  'DGVsV1';
  'HpVsV1';
  'CAVsVIS';
  'DGVsVIS';
  'HpVsVIS';
  'DGVsCA';
  'LPVsLGN';
  'LPVsV1';
  'LPVsV2';
  'LPVsVIS';
  'LPVsCA';
  'LPVsDG';
  'LPVsHp';
  'LGNVsLP';
  'V1VsLP';
  'V2VsLP';
  'VISVsLP';
  'CAVsLP';
  'DGVsLP';
  'HpVsLP';
  'LGNVsV2';
  'ThVsV2';
  'V1VsV2';
  'V2VsV1';
  'V2VsVIS';
  'V2VsCA';
  'V2VsDG';
  'V2VsHp';
  'V2VsLGN';
  'V2VsTh';
  'VISVsV2';
  'CAVsV2';
  'DGVsV2';
  'HpVsV2';
  'PoVsVB';
  'PoVsS1';
  'PoVsRSC';
  'PoVsCA';
  'PoVsDG';
  'PoVsHp'
  'VBVsPo';
  'S1VsPo';
  'RSCVsPo';
  'CAVsPo';
  'DGVsPo';
  'HpVsPo';
  'CxVsTh';
  'CxVsHp';
  'ThVsCx';
  'HpVsCx'};

iAreas2compareOI = find(ismember(areas2compare,areas2compareCritical));

areas2compareCriticalReduced = {
  'ThVsS1';
  'ThVsRSC';
  'ThVsV1';
  'ThVsVIS';
  'ThVsCA';
  'ThVsDG';
  'ThVsHp';
  'S1VsRSC';
  'S1VsCA';
  'S1VsDG';
  'S1VsHp';
  'RSCVsCA';
  'RSCVsDG';
  'RSCVsHp';
  'V1VsCA';
  'V1VsDG';
  'V1VsHp';
  'VISVsCA';
  'VISVsDG';
  'VISVsHp';
  'CAVsDG';
  'S1VsTh';
  'RSCVsTh';
  'V1VsTh';
  'VISVsTh';
  'CAVsTh';
  'DGVsTh';
  'HpVsTh';
  'RSCVsS1';
  'CAVsS1';
  'DGVsS1';
  'HpVsS1';
  'CAVsRSC';
  'DGVsRSC';
  'HpVsRSC';
  'CAVsV1';
  'DGVsV1';
  'HpVsV1';
  'CAVsVIS';
  'DGVsVIS';
  'HpVsVIS';
  'DGVsCA';
  'ThVsV2';
  'V1VsV2';
  'V2VsV1';
  'V2VsCA';
  'V2VsDG';
  'V2VsHp';
  'V2VsTh';
  'CAVsV2';
  'DGVsV2';
  'HpVsV2';
  'CxVsTh';
  'CxVsHp';
  'ThVsCx';
  'HpVsCx'};

iAreas2compareOIReduced = find(ismember(areas2compare,areas2compareCriticalReduced));

areasCohMat = {
  'VB';
  'LGN';
  'Th';
  'lS1';
  'V1';
  'lRSC';
  'DG';
  'CA1';
  'CA3';
  'CA';
  'Hp'};

awake = {
  '20181005134915';
  '20181008164836';
  '20181019140918';
  '20181008184928';
  '20181128182010';
  '20181129171747';
  '20181129195448';
  '20181219182543';
  '20181221145727';
  '20181221165147';
  '20190102134858';
  '20190102162201';
  '20190106174515';
  '20190122191103';
  '20190205085751';
  '20190208100404';
  '20190212095950';
  '20190215112507';
  '20190222104514';
  '20190312165540';
  '20190208152417';
  '20190211091018';
  '20190218132050';
  '20190311150750';
  '20190214143752';
  '20190218173940';
  '20190221114741';
  '20190225102939';
  '20190312133144';
  '20190315110943';
  '20190326130324';
  '20190331160303';
  '20190404111519';
  '20190405102743';
  '20190511141202';
  '20190514115945';
  '20190514140920';
  '20190515122749';
  '20190515124626';
  '20190515135343';
  '20190521102142';
  '20190522131946';
  '20190601133104';
  '20190602202843';
  '20190604174736';
  '20190606142912';
  '20190617184125';
  '20190621114337';
  '20190601100927';
  '20190602180606';
  '20190604103515';
  '20190606105826';
  '20190618152633';
  '20190621174558';
  '20190624192349';
  '20190624202957';
  '20191202210205';
  '20191203203543';
  '20191204184925';
  '20191205204947';
  '20191207175831';
  '20191208173908';
  '20191204233938';
  '20191206125358';
  '20191207205552';
  '20191209165319';
  '20191212205436';
  '20191213153806';
  '20191214153637';
  '20191212233211';
  '20191213202419';
  '20191215192006';
  '20191216155448';
  '20191217162628';
  '00000766640955';
  '00000767871931';
  '00000768515987';
  '00000771160300';
  '00000771990200';
  '00000774875821';
  '00000778240327';
  '00000778998620';
  '00000779839471';
  '00000781842082';
  '00000786091066';
  '00000787025148';
  '00000789848216';
  '00000793224716';
  '00000794812542';
  '00000816200189';
  '00000819186360';
  '00000819701982';
  '00000821695405';
  '00000829720705';
  '00000831882777';
  '00000835479236';
  '00000839068429';
  '00000839557629';
  '00000840012044';
  '00000847657808'};

awakePaired = {
  '20181128182010';
  '20181129171747';
  '20181221145727';
  '20190102134858';
  '20190106174515';
  '20190208100404';
  '20190212095950';
  '20190215112507';
  '20190211091018';
  '20190218132050';
  '20190311150750';
  '20190218173940';
  '20190221114741';
  '20190225102939';
  '20190514115945';
  '20190522131946';
  '20190604174736';
  '20190606142912';
  '20190617184125';
  '20190604103515';
  '20190606105826';
  '20190618152633';
  '20191203203543';
  '20191208173908';
  '20191209165319';
  '20191214153637';
  '20191217162628'};

anaesthesia = {
  '20180313172457';
  '20180711';
  '20180712';
  '20180713';
  '20180717';
  '20180719';
  '20181128204254';
  '20181129195448';
  '20181221165147';
  '20190102162201';
  '20190106200423';
  '20181128204254';
  '20181129195448';
  '20181221165147';
  '20190102162201';
  '20190106200423';
  '20190208110536';
  '20190212111502';
  '20190215123537';
  '20190211102014';
  '20190218142753';
  '20190311161859';
  '20190218185253';
  '20190221125927';
  '20190225114909';
  '20190514140920';
  '20190522153408';
  '20190604200258';
  '20190606163816';
  '20190617204956';
  '20190604124633';
  '20190606124357';
  '20190618173028';
  '20191018163127';
  '20191106163201';
  '20191107210153';
  '20191203224555';
  '20191208195945';
  '20191209191616';
  '20191214174709';
  '20191217200806';
  '20200316232145';
  '20200317192245';
  '20200318174631';
  '20200319130824';
  '20200323122038';
  '20200324161349';
  '20200325122741';
  '20210802135753'};

anaesthesiaPaired = {
  '20181128204254';
  '20181129195448';
  '20181221165147';
  '20190102162201';
  '20190106200423';
  '20190208110536';
  '20190212111502';
  '20190215123537';
  '20190211102014';
  '20190218142753';
  '20190311161859';
  '20190218185253';
  '20190221125927';
  '20190225114909';
  '20190514140920';
  '20190522153408';
  '20190604200258';
  '20190606163816';
  '20190617204956';
  '20190604124633';
  '20190606124357';
  '20190618173028';
  '20191203224555';
  '20191208195945';
  '20191209191616';
  '20191214174709';
  '20191217200806'};

% List of series to discard
except = {
  '201902050857512'; % Probe missed thalamus of M190128_A_MD (this and below)
  '201902050857513';
  '201902050857514';
  '2019020508575124';
  '201902081004042';
  '201902081004043';
  '201902081004044';
  '2019020810040424';
  '201902081105362';
  '201902081105363';
  '201902081105364';
  '2019020811053624';
  '201902120959502';
  '201902120959503';
  '201902120959504';
  '2019021209595024';
  '201902121115022';
  '201902121115023';
  '201902121115024';
  '2019021211150224';
  '201902151125072';
  '201902151125073';
  '201902151125074';
  '2019021511250724';
  '201902151235372';
  '201902151235373';
  '201902151235374';
  '2019021512353724';
  '201902221045142';
  '201902221045143';
  '201902221045144';
  '2019022210451424';
  '201903121655402';
  '201903121655403';
  '201903121655404';
  '2019031216554024';
  '201906242029571'; % extra recording in M190523_B_MD using Todor's ePhys box (this and below)
  '201906242029572';
  '201906242029573';
  '201906242029574';
  '2019062420295724';
  '201906242029575';
  '201906242029576';
  '2019062420295756';
  '201906242029577';
  '201905211021422'; % Probe missed VB in M190503_B_MD (this and below)
  '201905221319462';
  '201905221534082';
  '201912032035431'; % LFP trace saturations in multiple Neuropixels recordings (this and below)
  '201912032035432';
  '201912032035433';
  '201912032035434';
  '2019120320354324';
  '201912032035435';
  '201912032035436';
  '2019120320354356';
  '201912032035437';
  '201912032035438';
  '201912032035439';
  '2019120320354310';
  '20191203203543810';
  '2019120320354311';
  '201912041849251';
  '201912041849252';
  '201912041849253';
  '201912041849254';
  '2019120418492524';
  '201912041849255';
  '201912041849256';
  '2019120418492556';
  '201912041849257';
  '201912041849258';
  '201912041849259';
  '2019120418492510';
  '20191204184925810';
  '2019120418492511';
  '201912052049472';
  '201912052049473';
  '201912052049474';
  '2019120520494724';
  '201912052049475';
  '201912052049476';
  '2019120520494756';
  '201912052049477';
  '201912141536372';
  '201912141536373';
  '201912141536374';
  '2019121415363724';
  '201912141536375';
  '201912141536376';
  '2019121415363756';
  '201912141536377';
  '201912141747092';
  '201912141747093';
  '201912141747094';
  '2019121417470924';
  '201912141747095';
  '201912141747096';
  '2019121417470956';
  '201912141747097';
  '201912122332111';
  '201912122332118';
  '201912122332119';
  '2019121223321110';
  '20191212233211810';
  '2019121223321111';
  '201912161554482';
  '201912161554483';
  '201912161554484';
  '2019121615544824';
  '201912161554485';
  '201912161554486';
  '2019121615544856';
  '201912161554487'};

exceptFR = {
  '201902050857512'; % Probe missed thalamus of M190128_A_MD (this and below)
  '201902050857513';
  '201902050857514';
  '2019020508575124';
  '201902081004042';
  '201902081004043';
  '201902081004044';
  '2019020810040424';
  '201902081105362';
  '201902081105363';
  '201902081105364';
  '2019020811053624';
  '201902120959502';
  '201902120959503';
  '201902120959504';
  '2019021209595024';
  '201902121115022';
  '201902121115023';
  '201902121115024';
  '2019021211150224';
  '201902151125072';
  '201902151125073';
  '201902151125074';
  '2019021511250724';
  '201902151235372';
  '201902151235373';
  '201902151235374';
  '2019021512353724';
  '201902221045142';
  '201902221045143';
  '201902221045144';
  '2019022210451424';
  '201903121655402';
  '201903121655403';
  '201903121655404';
  '2019031216554024';
  '201906242029571'; % extra recording in M190523_B_MD using Todor's ePhys box (this and below)
  '201906242029572';
  '201906242029573';
  '201906242029574';
  '2019062420295724';
  '201906242029575';
  '201906242029576';
  '2019062420295756';
  '201906242029577';
  '201905211021422'; % Probe missed VB in M190503_B_MD (this and below)
  '201905221319462';
  '201905221534082'};

except = exceptFR; % Include saturations

% except = {
%   '201902050857512'; % Probe missed thalamus of M190128_A_MD (this and below)
%   '201902050857513';
%   '201902050857514';
%   '2019020508575124';
%   '201902081004042';
%   '201902081004043';
%   '201902081004044';
%   '2019020810040424';
%   '201902081105362';
%   '201902081105363';
%   '201902081105364';
%   '2019020811053624';
%   '201902120959502';
%   '201902120959503';
%   '201902120959504';
%   '2019021209595024';
%   '201902121115022';
%   '201902121115023';
%   '201902121115024';
%   '2019021211150224';
%   '201902151125072';
%   '201902151125073';
%   '201902151125074';
%   '2019021511250724';
%   '201902151235372';
%   '201902151235373';
%   '201902151235374';
%   '2019021512353724';
%   '201902221045142';
%   '201902221045143';
%   '201902221045144';
%   '2019022210451424';
%   '201903121655402';
%   '201903121655403';
%   '201903121655404';
%   '2019031216554024';
%   '201906242029571'; % extra recording in M190523_B_MD using Todor's ePhys box (this and below)
%   '201906242029572';
%   '201906242029573';
%   '201906242029574';
%   '2019062420295724';
%   '201906242029575';
%   '201906242029576';
%   '2019062420295756';
%   '201906242029577';
%   '201905211021422'; % Probe missed VB in M190503_B_MD (this and below)
%   '201905221319462';
%   '201905221534082';
%   '201911072101532'; % LFP trace saturations in multiple Neuropixels recordings (this and below)
%   '201911072101533';
%   '201911072101534';
%   '2019110721015324';
%   '201911072101535';
%   '201911072101536';
%   '2019110721015356';
%   '201911072101537';
%   '201912032035431';
%   '201912032035432';
%   '201912032035433';
%   '201912032035434';
%   '2019120320354324';
%   '201912032035435';
%   '201912032035436';
%   '2019120320354356';
%   '201912032035437';
%   '201912032035438';
%   '201912032035439';
%   '2019120320354310';
%   '20191203203543810';
%   '2019120320354311';
%   '201912041849251';
%   '201912041849252';
%   '201912041849253';
%   '201912041849254';
%   '2019120418492524';
%   '201912041849255';
%   '201912041849256';
%   '2019120418492556';
%   '201912041849257';
%   '201912041849258';
%   '201912041849259';
%   '2019120418492510';
%   '20191204184925810';
%   '2019120418492511';
%   '201912052049472';
%   '201912052049473';
%   '201912052049474';
%   '2019120520494724';
%   '201912052049475';
%   '201912052049476';
%   '2019120520494756';
%   '201912052049477';
%   '201912072055521';
%   '201912072055528';
%   '201912072055529';
%   '2019120720555210';
%   '20191207205552810';
%   '2019120720555211';
%   '201912091653192';
%   '201912091653193';
%   '201912091653194';
%   '2019120916531924';
%   '201912091653195';
%   '201912091653196';
%   '2019120916531956';
%   '201912091653197';
%   '201912122054361';
%   '201912122054368';
%   '201912122054369';
%   '2019121220543610';
%   '20191212205436810';
%   '2019121220543611';
%   '201912141536372';
%   '201912141536373';
%   '201912141536374';
%   '2019121415363724';
%   '201912141536375';
%   '201912141536376';
%   '2019121415363756';
%   '201912141536377';
%   '201912141747092';
%   '201912141747093';
%   '201912141747094';
%   '2019121417470924';
%   '201912141747095';
%   '201912141747096';
%   '2019121417470956';
%   '201912141747097';
%   '201912122332111';
%   '201912122332118';
%   '201912122332119';
%   '2019121223321110';
%   '20191212233211810';
%   '2019121223321111';
%   '201912151920061';
%   '201912151920067';
%   '201912151920068';
%   '201912151920069';
%   '2019121519200610';
%   '20191215192006810';
%   '2019121519200611';
%   '201912161554482';
%   '201912161554483';
%   '201912161554484';
%   '2019121615544824';
%   '201912161554485';
%   '201912161554486';
%   '2019121615544856';
%   '201912161554487';
%   '201912171626281';
%   '201912171626288';
%   '201912171626289';
%   '2019121716262810';
%   '20191217162628810';
%   '2019121716262811';
%   '202003191308241';
%   '202003191308242';
%   '202003191308243';
%   '202003191308244';
%   '2020031913082424';
%   '202003191308245';
%   '202003191308246';
%   '2020031913082456';
%   '202003191308247';
%   '202003191308248';
%   '202003191308249';
%   '2020031913082410';
%   '20200319130824810';
%   '2020031913082411';
%   '202003241613491';
%   '202003241613498';
%   '202003241613499';
%   '2020032416134910';
%   '20200324161349810';
%   '2020032416134911'};

% Series containing sharp wave ripples
fullSeries = {'M190114_A_MD_s201901221911036';
                'M190128_A_MD_s201902050857516';
                'M190128_A_MD_s201902081004046';
                'M190128_A_MD_s201902120959506';
                'M190128_A_MD_s201902151125076';
                  'M190128_B_MD_s201902081524176';
                  'M190128_B_MD_s201902181320506';
                  'M190128_B_MD_s201903111507506';
                    'M190128_C_MD_s201902141437526';
                    'M190128_C_MD_s201903121331446';
                      'M190322_C_MD_s201904041115196';
                      'M190322_C_MD_s201904051027436';
                        'M190503_B_MD_s201905221319466';
                          'M190523_A_MD_s201906011331046';
                          'M190523_A_MD_s201906022028436';
                          'M190523_A_MD_s201906041747366';
                          'M190523_A_MD_s201906061429126';
                          'M190523_A_MD_s201906171841256';
                          'M190523_A_MD_s201906211143376';
                            'M190523_B_MD_s201906021806066';
                            'M190523_B_MD_s201906041035156';
                            'M190523_B_MD_s201906061058266';
                            'M190523_B_MD_s201906181526336';
                            'M190523_B_MD_s201906211745586';
                              'M191018_MD_s201910181631276';
                                'M191106_MD_s201911061632016';
                                  'M191119_A_MD_s201912022102056';
                                  'M191119_A_MD_s201912032245556';
                                  	'M191119_B_MD_s2019120717583111';
                                    'M191119_B_MD_s2019120817390811';
                                    	'M191119_C_MD_s2019120423393811';
                                      'M191119_C_MD_s201912061253586';
                                      	'M191128_A_MD_s201912131538066';
                                        	'M191128_B_MD_s201912132024196';
                                          'M191128_B_MD_s2019121720080611';
                                          	'M200316_MD_s202003162321456';
                                            	'M200317_MD_s202003171922456';
                                              	'M200318_MD_s202003181746316';
                                                	'M200323_MD_s202003231220386';
                                                  	'M200325_MD_s2020032512274111'};

% Order of channels with sharp wave ripples in fullSeries
ch = {1;          % M190114_A_MD
      1;1;1;1;    % M190128_A_MD
      1;1;1;      % M190128_B_MD
      1;1;        % M190128_C_MD
      1;1;        % M190322_C_MD
      1;          % M190503_B_MD
      2;2;2;2;2;2;% M190523_A_MD
      1;1;1;1;1;  % M190523_B_MD
      1;          % M191018_MD
      1;          % M191106_MD
      2;1;        % M191119_A_MD
      3;2;        % M191119_B_MD
      3;2;        % M191119_C_MD
      1;          % M191128_A_MD
      1;2;        % M191128_B_MD
      1;          % M200316_MD
      1;          % M200317_MD
      1;          % M200318_MD
      1;          % M200323_MD
      2};         % M200325_MD

% Positive and negative subpopulation peaks and troughs for each brain
% area for Allen data. Each cell is indexing into phaseHistoBinCentres
% variable defined in params.m. The first value of each array is the
% location of the positive peak, the second value is the location of a
% trough to the right, the third value is the location of the negative
% peak, and the fourth value is the location of the trough to the right.
% The peak and trough locations are derived on MUA histograms (not unit
% histograms) at 0.3 Hz.
modesAllensdk = {
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [04 07 10 16]; % rV1
  [04 08 12 16]; %[04 10 12 15]; % rVB
  [];
  [04 08 12 16]; %[04 11 13 15]; % rLP
  [04 08 12 16]; %[04 12 15 16]; % rLGN
  [04 08 12 16]; %[03 11 15 16]; % rTh
  [06 07 10 01]; % rDG
  [05 08 11 01]; % rCA1
  [04 08 11 16]; % rCA3
  [05 08 11 16]; % rCA
  [05 08 11 16]; % rHp
  [04 08 12 16]; %[04 10 12 15]; % VB
  [04 08 12 16]; %[04 12 15 16]; % LGN
  [04 08 12 16]; %[03 11 15 16]; % Th
  [06 07 10 01]; % DG
  [05 08 11 01]; % CA1
  [04 08 11 16]; % CA3
  [05 08 11 16]; % CA
  [05 08 11 16]; % Hp
  [04 07 10 16]; % V1
  [];
  [];
  [04 07 10 15]; % rVIS
  [04 07 10 15]; % VIS
  [];
  [];
  [04 08 12 16]; %[04 11 13 15]; % LP
  [];
  [04 07 10 15]; % rV2
  [04 07 10 15]; % V2
  []};

% Positive and negative subpopulation peaks and troughs for each brain
% area for UOL data. The peak and trough locations are derived on MUA
% histograms at 0.03 Hz.
modesUOL = {
  [01 06 09 13]; % lS1
  [01 05 09 13]; % lVB1
  [01 05 09 13]; % lPo
  [01 05 09 13]; % lLP1
  [01 05 09 13]; % lTh1
  [01 05 09 14]; % lDG
  [02 05 09 14]; % lCA1
  [01 05 09 14]; % lHp1
  [01 05 09 13]; % lRSC
  [01 05 09 13]; % lVB2
  [01 05 09 13]; % lLP2
  [01 05 09 13]; % lLGN
  [01 05 09 12]; % lTh2
  [01 05 09 15]; % lCA3
  [01 05 09 12]; % lVB
  [01 05 09 13]; % lLP
  [01 05 09 11]; % lTh
  [01 05 09 14]; % lCA
  [01 05 09 14]; % lHp
  [02 04 07 13]; % lmPFC
  [02 04 08 13]; % lV1
  [02 05 08 13]; % rV1
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [];
  [01 05 09 12]; % VB
  [01 05 09 13]; % LGN
  [01 05 09 11]; % Th
  [01 05 09 14]; % DG
  [02 05 09 14]; % CA1
  [01 05 09 15]; % CA3
  [01 05 09 14]; % CA
  [01 05 09 14]; % Hp
  [];
  [01 06 09 13]; % Cx
  [];
  [];
  [];
  [01 06 09 13]; % S1
  [01 05 09 13]; % RSC
  [01 05 09 13]; % LP
  [];
  [];
  [];
  [01 05 09 13]}; % Po