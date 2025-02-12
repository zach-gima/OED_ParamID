%% Reference Potential for Neg Electrode: Unref(theta_n)
%   Created May 5, 2017 by Saehong Park & Dylan Kato
%   NCA cell, (NCR-18650pf) 
%   Curve-fitting w/ Bosch OCP data

function [Uref,varargout] = refPotentialAnode(p,theta)


 c_n=[-0.084294960339275
   0.920754744005144
  -0.500066623566425
   0.062731837918546
   0.782151587417570
  -0.373761901864611
   0.019988184317997
   0.543282314780430
  -0.295609630222051
   0.040970248093866
   0.231152288743602
  -0.217847875913234
   0.068744203951316
   0.353848415118256
  -0.114753994434564
  -0.028613032233089
   0.260671608316041
  -0.212058177468640
  -0.025506157489854
   0.211741908826122
  -0.241880220004548
   0.188872027034948
   0.867520021192469
  -0.225038983698359
  -0.111904175370177
   0.537399173641857
  -0.020780743382893
   0.108353745941168
   0.537735904911254
  -0.020226723056513
   0.171375773597772
   0.729717193288193
  -0.323902793428930
   0.064143152162965
   1.289849595601526
   0.704961322280748
   0.023028206444624
   0.481699223765299
  -0.076233450161839
  -0.182559256738691
   0.830851470359638
  -0.226362977193547
  -0.040952011143767
   1.626936110900125
   0.295695270567609
  -1.000228763094078
   0.007914258576845
  -0.016476666187381
  -0.341740372496750
   0.001274961492701
  -0.004879090290810
  -0.930906698538900
   0.001549868904555
  -0.010583717929547
   2.554274538083029
  -0.012402969675540
  -0.029257893810540
  -0.512533408582419
   0.066122834568301
  -0.077930639597751
  -0.499673574757569
   0.044470609922510
  -0.134483437256594
   1.904111886758372
  -0.035336812622768
  -0.306171040837701
  -1.122974595772499
   0.028740372472439
  -0.079271479637875
  -0.093855421675871
   0.930843806570863
  -0.516652668839875
  -0.846383609865041
   0.012151749801329
  -0.029511731110250
  -0.561782895480513
   0.098392530745244
  -0.109853910868333
  -0.818206413176353
   0.026850808833446
  -0.051805538572186
  -0.525543070925015
   0.188590232596615
  -0.192054642003214
  -0.046580230674248
   0.002863828671823
  -0.000914487593373
   2.650656293235332
  -0.008182255230700
  -0.117937922743741
  -0.295664205008775
   0.137690106957231
  -0.310460986123659
  -0.835065551163236
   0.711574616090746
  -0.997353098073145
   0.415746756470558
   0.423984781966332
   3.189835673119072
   0.413779708001205
   0.426343693564050
   3.190867502582611];
 
 Uref=c_n(1)*exp(-((theta - c_n(2)).^2/c_n(3)^2))+...
    c_n(4).*exp(-((theta - c_n(5)).^2/c_n(6).^2))+...
    c_n(7).*exp(-((theta - c_n(8)).^2/c_n(9).^2))+...
    c_n(10).*exp(-((theta - c_n(11)).^2/c_n(12).^2))+...
    c_n(13).*exp(-((theta - c_n(14)).^2/c_n(15).^2))+...
    c_n(16).*exp(-((theta - c_n(17)).^2/c_n(18).^2))+...
    c_n(19).*exp(-((theta - c_n(20)).^2/c_n(21).^2))+...
    c_n(22).*exp(-((theta - c_n(23)).^2/c_n(24).^2))+...
    c_n(25).*exp(-((theta - c_n(26)).^2/c_n(27).^2))+...
    c_n(28).*exp(-((theta - c_n(29)).^2/c_n(30).^2))+...
    c_n(31).*exp(-((theta - c_n(32)).^2/c_n(33).^2))+...
    c_n(34).*exp(-((theta - c_n(35)).^2/c_n(36).^2))+...
    c_n(37).*exp(-((theta - c_n(38)).^2/c_n(39).^2))+...
    c_n(40).*exp(-((theta - c_n(41)).^2/c_n(42).^2))+...
    c_n(43).*exp(-((theta - c_n(44)).^2/c_n(45).^2))+...
    c_n(46).*exp(-((theta - c_n(47)).^2/c_n(48).^2))+...
    c_n(49).*exp(-((theta - c_n(50)).^2/c_n(51).^2))+...
    c_n(52).*exp(-((theta - c_n(53)).^2/c_n(54).^2))+...
    c_n(55).*exp(-((theta - c_n(56)).^2/c_n(57).^2))+...
    c_n(58).*exp(-((theta - c_n(59)).^2/c_n(60).^2))+...
    c_n(61).*exp(-((theta - c_n(62)).^2/c_n(63).^2))+...
    c_n(64).*exp(-((theta - c_n(65)).^2/c_n(66).^2))+...
    c_n(67).*exp(-((theta - c_n(68)).^2/c_n(69).^2))+...
    c_n(70).*exp(-((theta - c_n(71)).^2/c_n(72).^2))+...
    c_n(73).*exp(-((theta - c_n(74)).^2/c_n(75).^2))+...
    c_n(76).*exp(-((theta - c_n(77)).^2/c_n(78).^2))+...
    c_n(79).*exp(-((theta - c_n(80)).^2/c_n(81).^2))+...
    c_n(82).*exp(-((theta - c_n(83)).^2/c_n(84).^2))+...
    c_n(85).*exp(-((theta - c_n(86)).^2/c_n(87).^2))+...
    c_n(88).*exp(-((theta - c_n(89)).^2/c_n(90).^2))+...
    c_n(91).*exp(-((theta - c_n(92)).^2/c_n(93).^2))+...
    c_n(94).*exp(-((theta - c_n(95)).^2/c_n(96).^2))+...
    c_n(97).*exp(-((theta - c_n(98)).^2/c_n(99).^2))+...
    c_n(100).*exp(-((theta - c_n(101)).^2/c_n(102).^2));
% Gradient of OCP wrt theta
if(nargout >= 2)

%     % Polynomial Fit
%     dUref = ppvalFast(p.dUppn,theta);
%     varargout{1} = dUref / p.c_s_n_max;

dUref=-2.*(theta - c_n(2))/(c_n(3).^2).*c_n(1).*exp(-((theta - c_n(2)).^2/c_n(3).^2))...
     -2.*(theta - c_n(5))/(c_n(6).^2).*c_n(4).*exp(-((theta - c_n(5)).^2/c_n(6).^2))...
     -2.*(theta - c_n(8))/(c_n(9).^2).*c_n(7).*exp(-((theta - c_n(8)).^2/c_n(9).^2))...
     -2.*(theta - c_n(11))/(c_n(12).^2).*c_n(10).*exp(-((theta - c_n(11)).^2/c_n(12).^2))...
     -2.*(theta - c_n(14))/(c_n(15).^2).*c_n(13).*exp(-((theta - c_n(14)).^2/c_n(15).^2))...
     -2.*(theta - c_n(17))/(c_n(18).^2).*c_n(16).*exp(-((theta - c_n(17)).^2/c_n(18).^2))...
     -2.*(theta - c_n(20))/(c_n(21).^2).*c_n(19).*exp(-((theta - c_n(20)).^2/c_n(21).^2))...
     -2.*(theta - c_n(23))/(c_n(24).^2).*c_n(22).*exp(-((theta - c_n(23)).^2/c_n(24).^2))...
     -2.*(theta - c_n(26))/(c_n(27).^2).*c_n(25).*exp(-((theta - c_n(26)).^2/c_n(27).^2))...
     -2.*(theta - c_n(29))/(c_n(30).^2).*c_n(28).*exp(-((theta - c_n(29)).^2/c_n(30).^2))...
     -2.*(theta - c_n(32))/(c_n(33).^2).*c_n(31).*exp(-((theta - c_n(32)).^2/c_n(33).^2))...
     -2.*(theta - c_n(35))/(c_n(36).^2).*c_n(34).*exp(-((theta - c_n(35)).^2/c_n(36).^2))...
     -2.*(theta - c_n(38))/(c_n(39).^2).*c_n(37).*exp(-((theta - c_n(38)).^2/c_n(39).^2))...
     -2.*(theta - c_n(41))/(c_n(42).^2).*c_n(40).*exp(-((theta - c_n(41)).^2/c_n(42).^2))...
     -2.*(theta - c_n(44))/(c_n(45).^2).*c_n(43).*exp(-((theta - c_n(44)).^2/c_n(45).^2))...
     -2.*(theta - c_n(47))/(c_n(48).^2).*c_n(46).*exp(-((theta - c_n(47)).^2/c_n(48).^2))...
     -2.*(theta - c_n(50))/(c_n(51).^2).*c_n(49).*exp(-((theta - c_n(50)).^2/c_n(51).^2))...
     -2.*(theta - c_n(53))/(c_n(54).^2).*c_n(52).*exp(-((theta - c_n(53)).^2/c_n(54).^2))...
     -2.*(theta - c_n(56))/(c_n(57).^2).*c_n(55).*exp(-((theta - c_n(56)).^2/c_n(57).^2))...
     -2.*(theta - c_n(59))/(c_n(60).^2).*c_n(58).*exp(-((theta - c_n(59)).^2/c_n(60).^2))...
     -2.*(theta - c_n(62))/(c_n(63).^2).*c_n(61).*exp(-((theta - c_n(62)).^2/c_n(63).^2))...
     -2.*(theta - c_n(65))/(c_n(66).^2).*c_n(64).*exp(-((theta - c_n(65)).^2/c_n(66).^2))...
     -2.*(theta - c_n(68))/(c_n(69).^2).*c_n(67).*exp(-((theta - c_n(68)).^2/c_n(69).^2))...
     -2.*(theta - c_n(71))/(c_n(72).^2).*c_n(70).*exp(-((theta - c_n(71)).^2/c_n(72).^2))...
     -2.*(theta - c_n(74))/(c_n(75).^2).*c_n(73).*exp(-((theta - c_n(74)).^2/c_n(75).^2))...
     -2.*(theta - c_n(77))/(c_n(78).^2).*c_n(76).*exp(-((theta - c_n(77)).^2/c_n(78).^2))...
     -2.*(theta - c_n(80))/(c_n(81).^2).*c_n(79).*exp(-((theta - c_n(80)).^2/c_n(81).^2))...
     -2.*(theta - c_n(83))/(c_n(84).^2).*c_n(82).*exp(-((theta - c_n(83)).^2/c_n(84).^2))...
     -2.*(theta - c_n(86))/(c_n(87).^2).*c_n(85).*exp(-((theta - c_n(86)).^2/c_n(87).^2))...
     -2.*(theta - c_n(89))/(c_n(90).^2).*c_n(88).*exp(-((theta - c_n(89)).^2/c_n(90).^2))...
     -2.*(theta - c_n(92))/(c_n(93).^2).*c_n(91).*exp(-((theta - c_n(92)).^2/c_n(93).^2))...
     -2.*(theta - c_n(95))/(c_n(96).^2).*c_n(94).*exp(-((theta - c_n(95)).^2/c_n(96).^2))...
     -2.*(theta - c_n(98))/(c_n(99).^2).*c_n(97).*exp(-((theta - c_n(98)).^2/c_n(99).^2))...
     -2.*(theta - c_n(101))/(c_n(102).^2).*c_n(100).*exp(-((theta - c_n(101)).^2/c_n(102).^2));
 
    varargout{1} = dUref;

end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end

