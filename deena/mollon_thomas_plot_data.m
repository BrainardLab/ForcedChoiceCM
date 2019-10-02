% This script lists Rayleigh matching loci for dichromats with varying cone
% peak spectral sensitivity and optical density, as measured by Thomas and 
% Mollon (2003). The script takes data points which were obtained from
% plots in the paper with WebPlotDigitizer, using these points to calculate
% a slope and intercept for each condition's line of best fit. Data arrays 
% are saved as .mat variables for lambda max fit parameters and optical 
% density fit parameters. 


%% Lambda max data 
lambdaMaxFitParams = zeros(4, 2); 

% Lambda max of 531nm 
data_531 = [0.030043235165559568, 0.3663244671063328;
0.05491585619392503, 0.35770220940057756;
0.07832482837339078, 0.3495606586227956;
0.1017376940988966, 0.34160541402303146;
0.12514958643789237, 0.33360359287876284;
0.14856537232292827, 0.32578807791251196;
0.17197629127541397, 0.31773968022373894;
0.19538915700091986, 0.30978443562397473;
0.2188029961129357, 0.30187576756871504;
0.24221683522495158, 0.2939670995134553;
0.26599423419845963, 0.2858256947436498;
0.2886717064156456, 0.27794002294468934;
0.3124544590149589, 0.2700547891696583;
0.33586440458093475, 0.26195981493638076;
0.3592782436929504, 0.25405114688112107;
0.3826881892589263, 0.24595617264784359;
0.4061010549844321, 0.2380009280480794;
0.42951100055040786, 0.22990595381480186;
0.45292483966242375, 0.22199728575954214;
0.47633770538792947, 0.21404204115977798;
0.49974667756739516, 0.20590049038199604;
0.5231605166794111, 0.1979918223267363;
0.5465704622453869, 0.18989684809345878;
0.5699852747439127, 0.18203475658270352;
0.5875510077037296, 0.17635942653603315;
0.9562892854294234, 0.05037732005830681;
0.979701177768419, 0.042375498914038245;
0.9958009906441463, 0.0370721971914974];
lambdaMaxFitParams(1,:) = polyfit(data_531(:,1), data_531(:,2), 1);

% Lambda max of 541nm 
data_541 = [0.03394787191254073, 0.2710961180580775;
0.06622147503875658, 0.26402933199533385;
0.09849366233004886, 0.25689479823149286;
0.13076726545626477, 0.24982801216874928;
0.16303803691263324, 0.24262573070381085;
0.1953095162864637, 0.23545732308942122;
0.2275824114952178, 0.2283566631761289;
0.25985318295158627, 0.22115438171119056;
0.29212725802277667, 0.21411017821547937;
0.3243961416992471, 0.20681756648241115;
0.35667186857784844, 0.1998524019713136;
0.3889405162818315, 0.19254849895472914;
0.42121057982073823, 0.18531234363924215;
0.45348559878187783, 0.17831330527759587;
0.858370024122499, 0.08953305782631554;
0.890644807111151, 0.08252272818115314;
0.9229141627325961, 0.07525269901511739;
0.9551898896111972, 0.06828753450401981;
0.9859966523606232, 0.061548740759692466]; 
lambdaMaxFitParams(2,:) = polyfit(data_541(:,1), data_541(:,2), 1);

% Lambda max of 551nm 
data_551 = [0.03268494716063142, 0.21066516867921736;
0.06502934203303051, 0.2069857676713426;
0.09737302898796785, 0.20327249281291912;
0.12971742386036705, 0.19959309180504436;
0.16206323456768984, 0.19598143849826694;
0.19440550568770348, 0.1922004159387461;
0.22675202431248814, 0.1886226364825174;
0.2590964191848873, 0.18494323547464264;
0.2914415219747482, 0.18129770831731656;
0.32378662476460923, 0.17765218115999043;
0.3502566194076557, 0.17494794656889662;
0.7354505396647135, 0.1312052917386785;
0.767797058289498, 0.12762751228244978;
0.8001393294095118, 0.12384648972292894;
0.8324872638692202, 0.12033645796779763;
0.8648330745765429, 0.1167248046610202;
0.8971781773664038, 0.11307927750369412;
0.9295246959911887, 0.1095014980474654;
0.9618697987810494, 0.10585597089013932;
0.988331500676686, 0.10275492833547512];
lambdaMaxFitParams(3,:) = polyfit(data_551(:,1), data_551(:,2), 1);

% Lambda max of 561nm 
data_561 = [0.03186659457475299, 0.17150699744493358;
0.06428815245049185, 0.17151984614686583;
0.0967111261611544, 0.17160044254989548;
0.12913551570674064, 0.17174878665402243;
0.16156132108725052, 0.17196487845924685;
0.19398641855029858, 0.17214709641392256;
0.22641647143557947, 0.1725664313224391;
0.24852437980987968, 0.1726852818153125;
0.4967451707361414, 0.1730897239102257;
0.5324735543451102, 0.17351197897827225;
0.5634265964393994, 0.1737726518856556;
0.590652495235084, 0.1736592384517167;
0.6272888925475082, 0.17387699891952968;
0.6597125741756324, 0.17399146917310804;
0.6921348399688334, 0.17403819172558893;
0.7245571057620339, 0.17408491427806994;
0.7569836190600057, 0.174334879933843;
0.7894073006881299, 0.1744493501874213;
0.8218309823162544, 0.1745638204409996;
0.8542589114491499, 0.17488153379787014;
0.8866804693248886, 0.17489438249980238;
0.9191020272006276, 0.17490723120173463;
0.9515242929938283, 0.17495395375421563;
0.9839515142092619, 0.1752377932605374];
lambdaMaxFitParams(4,:) = polyfit(data_561(:,1), data_561(:,2), 1);

save('mollon_thomas_lambdaMaxFitParams.mat', 'lambdaMaxFitParams'); 

%% Optical density data 
odFitParams = zeros(5, 2); 

% Optical density of 0.05, or -83.33%
data_5 = [0.031327438583960615, 0.15770117424892546
0.06375133973315256, 0.1567268263116374
0.096172328555894, 0.15570739070701778
0.12859622970508594, 0.15473304276972974
0.1610192987610063, 0.1537458126417755
0.18460717328768805, 0.15314333066663585
0.32018501207507283, 0.14886720471878945
0.3526022564780921, 0.14778979925617208
0.38502490948737667, 0.14679612803288478
0.41744714645002545, 0.14579601571426443
0.4498710475992175, 0.14482166777697636
0.48229328456186626, 0.14382155545835598
0.514715521524515, 0.1428214431397356
0.5471399774025547, 0.14185568332955834
0.5795599954498126, 0.14082121850249485
0.6119820243891435, 0.13981788563620795
0.6451449630457095, 0.13885437964311154
0.6768271223843948, 0.13782732264196682
0.709248527253772, 0.1368143281326803
0.7416699321231492, 0.13580133362339375
0.7740925851324337, 0.13480766240010644
0.8065135739551752, 0.13378822679548685
0.8389370590577314, 0.1328074377628657
0.8713601281136518, 0.13182020763491148
0.9037807008897574, 0.1307943309349588
0.9362021057591348, 0.12978133642567224
0.9686247587684194, 0.12878766520238494
0.9922046733713075, 0.1280619503403822]; 
odFitParams(1,:) = polyfit(data_5(:,1), data_5(:,2), 1);

% Optical density of 0.2, or -33.33% 
data_20 = [0.29491266136245664, 0.15082484483328285;
0.3261961925511709, 0.15048477392981813;
0.358659202084127, 0.15011588895383965;
0.3911234597569906, 0.14976632726386038;
0.42358896556976144, 0.14943608885988038;
0.4560536392892608, 0.14909296826523422;
0.4885170648688527, 0.1487305243845888;
0.5209817385883521, 0.14838740378994264;
0.5534459962612155, 0.1480378420999634;
0.5859098378874433, 0.14768183931465106;
0.6183745116069426, 0.14733871872000492;
0.6508404334663492, 0.147014921411358;
0.683304275092577, 0.14665891862604566;
0.7157681167188048, 0.14630291584073335;
0.7482336225315757, 0.14597267743675335;
0.7806970481111676, 0.14561023355610792;
0.813159641597488, 0.14523490748479637;
0.8456284757833452, 0.14495619784348102;
0.8780927334562088, 0.14460663615350175;
0.9105565750824365, 0.14425063336818944;
0.9430204167086642, 0.1438946305828771;
0.9754855064747994, 0.143557951083564;
0.996145966361468, 0.14335697305317235]; 
odFitParams(2,:) = polyfit(data_20(:,1), data_20(:,2), 1); 

% Optical density of 0.35, or 16.67% 
data_35 = [0.2576617143366741, 0.15251103415879189;
0.28794625167880716, 0.1527078400792983;
0.32045544238833595, 0.15305391668529175;
0.35296130472477844, 0.15334846452862055;
0.3854692472943999, 0.15367521784861476;
0.4179771898640213, 0.15400197116860898;
0.4504830522004637, 0.15429651901193778;
0.4829914108167209, 0.15462971342726506;
0.515500185479614, 0.15496934893792544;
0.5480075733203876, 0.15528751413080885;
0.5805131582924062, 0.15557776791058228;
0.6130198527221203, 0.15588519794457725;
0.6455286273850133, 0.15622483345523763;
0.6780369860012705, 0.1565580278705649;
0.7105440964776204, 0.15687189899989296;
0.7430512069539702, 0.157185770129221;
0.7755599816168632, 0.15752540563988138;
0.8080662599999414, 0.15782639457854325;
0.8405729544296555, 0.15813382461253822;
0.8730817290925486, 0.15847346012319857;
0.9055909198020775, 0.15881953672919202;
0.9380971981851556, 0.1591205256678539;
0.9706047247081413, 0.15944083789251504;
0.9942485440742868, 0.15970447527367343];
odFitParams(3,:) = polyfit(data_35(:,1), data_35(:,2), 1); 
 
% Optical density of 0.5, or 66.67% 
data_50 = [0.23780015066710347, 0.15364956672442218;
0.27035080735799927, 0.154637605831946;
0.30290035459119974, 0.15560846868524825;
0.3354528141508506, 0.15662441920588205;
0.3680036095239584, 0.15761460534518357;
0.4005552369903377, 0.15861767367515123;
0.4331076965499886, 0.15963362419578503;
0.4656593240163679, 0.1606366925257527;
0.4982109514827472, 0.16163976085572035;
0.5307634110423981, 0.16265571137635418;
0.5633172574241683, 0.16369313221476495;
0.5958683301616998, 0.1646876124176218;
0.628418293441536, 0.16566491636625716;
0.6609715850944584, 0.16669374907755713;
0.6935236286074735, 0.1677032585028579;
0.7260765042137602, 0.1687256501188248;
0.7586285477267753, 0.16973515954412552;
0.7911810072864263, 0.17075111006475935;
0.8237326347528056, 0.171754178394727;
0.8562834301259132, 0.1727443645340285;
0.8888354736389283, 0.17375387395932923;
0.9213871011053076, 0.1747569422892969;
0.9539378964784155, 0.1757471284285984;
0.9864903560380662, 0.1767630789492322];
odFitParams(4,:) = polyfit(data_50(:,1), data_50(:,2), 1); 

% Optical density of 0.65, or 116.67% 
data_65 = [0.030492849032566965, 0.14478033701076598;
0.21422819642318008, 0.15449852121046456;
0.24385062443658173, 0.15587546160579738;
0.2720052306734373, 0.15738945532884435;
0.3016387625092741, 0.1589383014018445;
0.33127183669381155, 0.16048006226997824;
0.36090491087834903, 0.16202182313811198;
0.39025454635814427, 0.1637477566186506;
0.3772061450609854, 0.16291690822794705;
0.4098000706020033, 0.16457482125011122;
0.4498064216884582, 0.1666825317668452;
0.47944178412949245, 0.1682597186593109;
0.5090757736166287, 0.16981564993717743;
0.538708847801166, 0.1713574108053112;
0.5683419219857035, 0.17289917167344493;
0.5979763691241391, 0.17446218815617784;
0.6276094433086765, 0.1760039490243116;
0.6572425174932139, 0.17754570989244536;
0.6853946320729952, 0.17902112861121977;
0.713546136451044, 0.18048710039017235;
0.7431792106355815, 0.18202886125830608;
0.7728127424714182, 0.18357770733130624;
0.8024458166559557, 0.18511946819944;
0.832077517886595, 0.1866399734529746;
0.8617083038146357, 0.1881463082967764;
0.8913441239069693, 0.18973058039410848;
0.9209767404402074, 0.19126525605737585;
0.9506102722760443, 0.19281410213037598;
0.980245177065779, 0.1943842038179753;
0.9995036013945009, 0.19533875942290965];
odFitParams(5,:) = polyfit(data_65(:,1), data_65(:,2), 1); 

save('mollon_thomas_odFitParams.mat', 'odFitParams'); 