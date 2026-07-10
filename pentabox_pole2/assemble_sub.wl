(* Assemble the SUB-SECTOR contribution to poles(I[1]):  Sigma_{j != top} c_j(d) M_j.
   If [1/eps^2] == c_-2 = -J0/2, the (untabulated) top masters net-contribute ZERO -> poles are
   isolated to the sub-sectors (route closes; class-D finite -> from KNOWN masters). If not, the top
   masters DO contribute (route blocked by them); the mismatch quantifies their pole contribution.
   Masters: exact bubble-chains (bubblechain_laurents.wl, keyed by sector) + pySecDec
   (masterLaurents.wl, keyed by index) ; the 4 top (sector 255) masters are SKIPPED. *)
Off[General::stop];
red=Get["kira_target.m"]; rhs=red[[1,2]];
bcl=Get["analysis/bubblechain_laurents.wl"];      (* sector -> {m2,m1,m0} *)
psd=Get["masterLaurents.wl"];                     (* {idx} -> {m2,m1,m0,mp1} *)
sectorOf[idx_]:=Total[Table[If[idx[[i]]>0,2^(i-1),0],{i,8}]];
lau[idx_]:=Module[{s=sectorOf[idx],key=idx},
  If[s==255, Return["TOP"]];                       (* skip top masters *)
  If[KeyExistsQ[bcl,s], Return[bcl[s][[1]]/eps^2+bcl[s][[2]]/eps+bcl[s][[3]]]];
  If[KeyExistsQ[psd,key], With[{v=psd[key]}, Return[v[[1]]/eps^2+v[[2]]/eps+v[[3]]+v[[4]] eps]]];
  Missing[idx]];
terms=List@@rhs; miss={}; topcnt=0; Isub=0;
Do[
  pb=Cases[t,_pentabox,{0,Infinity}][[1]]; idx=List@@pb; cf=t/pb;
  L=lau[idx];
  Which[
   L==="TOP", topcnt++,
   Head[L]===Missing, AppendTo[miss,idx],
   True, Isub+=Series[(cf/.d->4-2eps)*L,{eps,0,1}]//Normal
  ];
,{t,terms}];
Print["top masters skipped: ",topcnt,"   masters missing (not yet evaluated): ",Length[miss]];
If[Length[miss]>0, Print["  first few missing: ",Take[miss,Min[6,Length[miss]]]]];
Is=Series[Isub,{eps,0,1}];
P3=SeriesCoefficient[Is,-3]; P2=SeriesCoefficient[Is,-2]; P1=SeriesCoefficient[Is,-1];
Print["[1/eps^3] sub-sector sum (expect ~0): ",N[P3,10]];
Print["[1/eps^2] = P2_sub = ",N[P2,12]];
Print["[1/eps]   = P1_sub = ",N[P1,12]];
cm2=-1.914807547604452136*^-5; cm1=2.197205645392600987*^-4;
Print["c_-2 = ",cm2,"   P2_sub - c_-2 = ",N[P2-cm2,8],"   (|.|/|c_-2| = ",N[Abs[P2-cm2]/Abs[cm2],6],")"];
Print["  => if P2_sub ~ c_-2 : poles isolated to sub-sectors (top net 0).  if O(1): top masters contribute."];
