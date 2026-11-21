const css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const NS="http://www.w3.org/2000/svg";
const el=(t,a={},p)=>{const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);p&&p.appendChild(e);return e};
const tip=document.getElementById("tip");
let R=[],E=null,F=[];
const f2=v=>v==null?"–":(+v).toFixed(2);

function showTip(ev,html){
  const t=ev.touches?ev.touches[0]:ev;tip.innerHTML=html;tip.hidden=false;
  let left=t.clientX+14+scrollX;if(t.clientX+tip.offsetWidth+30>innerWidth)left=t.clientX-tip.offsetWidth-14+scrollX;
  tip.style.left=left+"px";tip.style.top=(t.clientY+scrollY+14)+"px";
}
function hideTip(){tip.hidden=true}
function svgFor(box,W,H,label){box.innerHTML="";return el("svg",{viewBox:`0 0 ${W} ${H}`,role:"img","aria-label":label},box)}
function grid(svg,ticks,y,L,W,Rr,fmt){ticks.forEach(v=>{el("line",{x1:L,x2:W-Rr,y1:y(v),y2:y(v),stroke:css("--rule"),"stroke-width":1},svg);el("text",{x:L-8,y:y(v)+4,"text-anchor":"end"},svg).textContent=fmt?fmt(v):v})}
/* Dot with a bigger invisible hit target */
function dot(svg,cx,cy,r,fill,stroke,html,ring){
  const c=el("circle",{cx,cy,r,fill:ring?css("--bg"):fill,stroke:ring?fill:stroke,"stroke-width":2},svg);
  const h=el("circle",{cx,cy,r:Math.max(r+5,9),fill:"transparent"},svg);
  h.addEventListener("pointermove",e=>showTip(e,html));h.addEventListener("pointerdown",e=>showTip(e,html));h.addEventListener("pointerleave",hideTip);
  return c;
}
/* Crosshair hover for per-residue line charts */
function crosshair(box,svg,x,inv,top,bottom,fmt,dots){
  const line=el("line",{y1:top,y2:bottom,stroke:css("--ink-3"),"stroke-width":1,"stroke-dasharray":"3 3",visibility:"hidden"},svg);
  const marks=dots.map(d=>el("circle",{r:4,fill:d.color,stroke:css("--bg"),"stroke-width":2,visibility:"hidden"},svg));
  const hit=el("rect",{x:0,y:top,width:svg.viewBox.baseVal.width,height:bottom-top,fill:"transparent"},svg);
  const move=ev=>{
    const pt=svg.createSVGPoint();const t=ev.touches?ev.touches[0]:ev;pt.x=t.clientX;pt.y=t.clientY;
    const p=pt.matrixTransform(svg.getScreenCTM().inverse());const r=inv(p.x);if(!r)return;
    line.setAttribute("x1",x(r.n));line.setAttribute("x2",x(r.n));line.setAttribute("visibility","visible");
    dots.forEach((d,i)=>{const v=d.get(r);if(v==null){marks[i].setAttribute("visibility","hidden");return}marks[i].setAttribute("cx",x(r.n));marks[i].setAttribute("cy",d.y(v));marks[i].setAttribute("visibility","visible")});
    showTip(ev,fmt(r));
  };
  hit.addEventListener("pointermove",move);hit.addEventListener("pointerdown",move);
  hit.addEventListener("pointerleave",()=>{hideTip();line.setAttribute("visibility","hidden");marks.forEach(m=>m.setAttribute("visibility","hidden"))});
}
function axisRes(svg,x,y0,W,lo,hi){[lo,50,100,150,200,250,hi].forEach(n=>{el("text",{x:x(n),y:y0+16,"text-anchor":"middle"},svg).textContent=n});el("text",{x:W/2,y:y0+32,"text-anchor":"middle"},svg).textContent="residue number"}

/* Fig 2: confidence along the chain */
function drawConf(){
  const box=document.getElementById("chartConf"),W=800,H=230,L=40,Rr=12,T=12,B=H-42;
  const svg=svgFor(box,W,H,"AlphaFold confidence per residue: low for residues 1 to 29, near 100 for the rest");
  const x=n=>L+(n-1)/289*(W-L-Rr),y=v=>B-(v/100)*(B-T);
  el("rect",{x:x(1),y:T,width:x(29.5)-x(1),height:B-T,fill:css("--band")},svg);
  el("text",{x:x(15),y:T+14,"text-anchor":"middle",class:"lbl"},svg).textContent="signal peptide";
  grid(svg,[0,50,70,90,100],y,L,W,Rr);
  el("text",{x:W-Rr,y:y(90)-6,"text-anchor":"end"},svg).textContent="very high (90+)";
  el("path",{d:R.map((r,i)=>(i?"L":"M")+x(r[0]).toFixed(1)+" "+y(r[2]).toFixed(1)).join(""),fill:"none",stroke:css("--s1"),"stroke-width":2,"stroke-linejoin":"round"},svg);
  axisRes(svg,x,B,W,1,290);
  const inv=px=>{const n=Math.max(1,Math.min(290,Math.round(1+(px-L)/(W-L-Rr)*289)));const r=R[n-1];return {n,r}};
  crosshair(box,svg,x,inv,T,B,o=>`Residue ${o.n} · ${o.r[1]}<br>Confidence ${o.r[2].toFixed(1)}`,[{color:css("--s1"),get:o=>o.r[2],y}]);
}

/* Fig 3: distribution of pairwise RMSDs */
function drawDist(){
  const box=document.getElementById("chartDist"),W=800,H=210,L=150,Rr=16,T=10,B=H-40,max=0.5;
  const svg=svgFor(box,W,H,"Strip plot: AlphaFold vs lab RMSDs fall within the range of lab vs lab RMSDs");
  const x=v=>L+(v/max)*(W-L-Rr);
  [0,0.1,0.2,0.3,0.4,0.5].forEach(v=>{el("line",{x1:x(v),x2:x(v),y1:T,y2:B,stroke:css("--rule"),"stroke-width":1},svg);el("text",{x:x(v),y:B+16,"text-anchor":"middle"},svg).textContent=v.toFixed(1)});
  el("text",{x:(L+W-Rr)/2,y:B+32,"text-anchor":"middle"},svg).textContent="RMSD between two structures (Å)";
  const rows=[{label:"Different crystals",vals:E.lab_lab_values,col:css("--s2"),op:.55,y:T+30,med:E.lab_lab.median},{label:"Same crystal",vals:E.same_crystal_values,col:css("--s2"),op:.3,y:T+80,med:E.same_crystal.median},{label:"AlphaFold vs lab",vals:E.ai_lab_values,col:css("--s1"),op:.7,y:T+130,med:E.ai_lab.median}];
  rows.forEach(rw=>{
    el("text",{x:L-14,y:rw.y+4,"text-anchor":"end",class:"lbl"},svg).textContent=rw.label;
    let seed=7;const rnd=()=>{seed=(seed*16807)%2147483647;return seed/2147483647};
    rw.vals.forEach(v=>el("circle",{cx:x(v),cy:rw.y+(rnd()-.5)*28,r:3,fill:rw.col,"fill-opacity":rw.op},svg));
    el("line",{x1:x(rw.med),x2:x(rw.med),y1:rw.y-18,y2:rw.y+18,stroke:css("--ink"),"stroke-width":2},svg);
    el("text",{x:x(rw.med)+6,y:rw.y-12},svg).textContent="median "+rw.med.toFixed(3);
    const hit=el("rect",{x:L,y:rw.y-20,width:W-L-Rr,height:40,fill:"transparent"},svg);
    const h=e=>showTip(e,`${rw.label}<br>${rw.vals.length} pairs<br>median ${rw.med.toFixed(3)} Å<br>range ${Math.min(...rw.vals).toFixed(3)}–${Math.max(...rw.vals).toFixed(3)} Å`);
    hit.addEventListener("pointermove",h);hit.addEventListener("pointerdown",h);hit.addEventListener("pointerleave",hideTip);
  });
}

/* Fig 4: scatter, lab wobble (B-factor z) vs AI uncertainty */
function drawScatter(){
  const box=document.getElementById("chartScatter"),W=800,H=360,L=56,Rr=16,T=12,B=H-46;
  const svg=svgFor(box,W,H,"Scatter plot: residues that wobble more in crystals have lower AlphaFold confidence");
  const res=E.residues,xs=res.map(r=>r.zb),ys=res.map(r=>100-r.plddt);
  const x0=Math.floor(Math.min(...xs)),x1=Math.ceil(Math.max(...xs)),y1=Math.ceil(Math.max(...ys));
  const x=v=>L+(v-x0)/(x1-x0)*(W-L-Rr),y=v=>B-(v/y1)*(B-T);
  for(let v=x0;v<=x1;v++){el("line",{x1:x(v),x2:x(v),y1:T,y2:B,stroke:css("--rule"),"stroke-width":1},svg);el("text",{x:x(v),y:B+16,"text-anchor":"middle"},svg).textContent=v}
  for(let v=0;v<=y1;v++){el("line",{x1:L,x2:W-Rr,y1:y(v),y2:y(v),stroke:css("--rule"),"stroke-width":1},svg);el("text",{x:L-8,y:y(v)+4,"text-anchor":"end"},svg).textContent=v}
  el("text",{x:(L+W-Rr)/2,y:B+34,"text-anchor":"middle"},svg).textContent="lab wobble: crystal B-factor (standard deviations from the chain's average)";
  el("text",{x:14,y:(T+B)/2,"text-anchor":"middle",transform:`rotate(-90 14 ${(T+B)/2})`},svg).textContent="AI uncertainty (100 − pLDDT)";
  res.forEach(r=>dot(svg,x(r.zb),y(100-r.plddt),3.5,css("--s1"),css("--bg"),`Residue ${r.pos} · ${r.aa}<br>B-factor z ${r.zb.toFixed(2)}<br>Confidence ${r.plddt.toFixed(1)}`));
  [73,277,187].forEach(p=>{const r=res.find(q=>q.pos===p);if(r)el("text",{x:x(r.zb)+8,y:y(100-r.plddt)-6,class:"lbl"},svg).textContent=r.aa[0]+r.aa.slice(1).toLowerCase()+p});
  const c=E.correlations.uncertainty_vs_bfactor;
  el("text",{x:W-Rr,y:T+14,"text-anchor":"end",class:"lbl"},svg).textContent=`ρ = ${c.rho.toFixed(2)} (95% CI ${c.lo.toFixed(2)}–${c.hi.toFixed(2)})`;
}

/* Fig 5: wobble and AI error along the chain */
function drawWobble(){
  const box=document.getElementById("chartWobble"),W=800,H=240,L=40,Rr=12,T=12,B=H-42,lo=34,hi=290,max=0.8;
  const svg=svgFor(box,W,H,"Lab wobble and AlphaFold error per residue; both peak at residues 73, 187 and 276 to 278");
  const res=E.residues,x=n=>L+(n-lo)/(hi-lo)*(W-L-Rr),y=v=>B-(v/max)*(B-T);
  grid(svg,[0,0.2,0.4,0.6,0.8],y,L,W,Rr,v=>v.toFixed(1));
  el("text",{x:L+4,y:T+2,"dominant-baseline":"hanging"},svg).textContent="Å";
  const path=(key,col)=>{let d="",prev=null;res.forEach(r=>{d+=(prev!=null&&r.pos===prev+1?"L":"M")+x(r.pos).toFixed(1)+" "+y(Math.min(r[key],max)).toFixed(1);prev=r.pos});el("path",{d,fill:"none",stroke:col,"stroke-width":2,"stroke-linejoin":"round"},svg)};
  path("rmsf",css("--s2"));path("af_dev",css("--s1"));
  [[73,"Asn73"],[193,"Ser193 (AI only)"],[277,"loop 276–278"]].forEach(([p,t])=>{const r=res.find(q=>q.pos===p);if(r)el("text",{x:x(p),y:y(Math.max(r.rmsf,r.af_dev))-8,"text-anchor":"middle",class:"lbl"},svg).textContent=t});
  [lo,50,100,150,200,250,hi].forEach(n=>{el("text",{x:x(n),y:B+16,"text-anchor":"middle"},svg).textContent=n});
  el("text",{x:W/2,y:B+32,"text-anchor":"middle"},svg).textContent="residue number (residues missing from any crystal are left out)";
  const byPos={};res.forEach(r=>byPos[r.pos]=r);
  const inv=px=>{let n=Math.round(lo+(px-L)/(W-L-Rr)*(hi-lo));for(let k=0;k<6;k++){if(byPos[n+k])return{n:n+k,r:byPos[n+k]};if(byPos[n-k])return{n:n-k,r:byPos[n-k]}}return null};
  crosshair(box,svg,x,inv,T,B,o=>`Residue ${o.n} · ${o.r.aa}<br><i style="background:${css("--s2")}"></i>Lab wobble ${o.r.rmsf.toFixed(2)} Å<br><i style="background:${css("--s1")}"></i>AI error ${o.r.af_dev.toFixed(2)} Å`,
    [{color:css("--s2"),get:o=>Math.min(o.r.rmsf,max),y},{color:css("--s1"),get:o=>Math.min(o.r.af_dev,max),y}]);
}

/* Fig 6: Trp185 chi2 angles */
function drawTrp(){
  const box=document.getElementById("chartTrp"),W=800,H=230,L=20,Rr=20,T=10,B=H-44;
  const svg=svgFor(box,W,H,"Trp185 ring angles: most near 90 degrees, a few at 21, 44, 79, 113 and minus 92");
  const x=v=>L+(v+180)/360*(W-L-Rr),mid=(T+B)/2;
  [-180,-90,0,90,180].forEach(v=>{el("line",{x1:x(v),x2:x(v),y1:T,y2:B,stroke:css("--rule"),"stroke-width":1},svg);el("text",{x:x(v),y:B+16,"text-anchor":"middle"},svg).textContent=v+"°"});
  el("text",{x:W/2,y:B+32,"text-anchor":"middle"},svg).textContent="Trp185 ring angle χ2";
  const lab=E.trp185.filter(t=>t.name!=="AlphaFold"),af=E.trp185.find(t=>t.name==="AlphaFold");
  const sorted=[...lab].sort((a,b)=>a.chi2-b.chi2);const stack={};
  /* Dots stack upward from the baseline in 3° bins; AlphaFold is a ring on the baseline row below */
  const base=B-30;
  sorted.forEach(t=>{const bin=Math.round(t.chi2/3);stack[bin]=(stack[bin]||0)+1;const k=stack[bin]-1;
    dot(svg,x(t.chi2),base-k*8,3.5,css("--s2"),css("--bg"),`${t.name}<br>χ2 ${t.chi2.toFixed(1)}°<br>released ${t.released}`)});
  if(af){dot(svg,x(af.chi2),B-10,5,css("--s1"),css("--bg"),`AlphaFold<br>χ2 ${af.chi2.toFixed(1)}°`,true);el("text",{x:x(af.chi2)+12,y:B-6},svg).textContent="AlphaFold"}
  const lbl=(name,dx,anchor)=>{const t=lab.find(q=>q.name===name);if(t)el("text",{x:x(t.chi2)+dx,y:base-14,"text-anchor":anchor,class:"lbl"},svg).textContent=name}
  lbl("5XG0-A",0,"middle");lbl("6ANE-C",0,"middle");lbl("6EQG-B",0,"middle");lbl("5XG0-C",10,"start");
}

/* Fig 7: blind test */
function drawBlind(){
  const fam=[...F].sort((a,b)=>a.first_structure.localeCompare(b.first_structure));
  const rows=[{label:"PETase (5XJH)",full:0.29,core:0.29,blind:false,year:"2018"}].concat(fam.map(f=>({label:`${f.organism.replace(/ DSM.*| TB100| YX/,"").replace("unidentified prokaryotic organism","compost metagenome (LCC)").replace(/^unidentified$/,"metagenome")} (${f.pdb.split("-")[0]})`,full:f.af_vs_lab_rmsd,core:f.af_core_rmsd,blind:f.blind,year:f.first_structure.slice(0,4),off:f.af_off_by_2A,offpl:f.af_off_mean_plddt,strict:f.strict_blind,mdate:f.af_model_date,pct:f.af_core_pct})));
  const rowH=26,T=10,B=T+rows.length*rowH+12,W=800,L=290,Rr=16,H=B+40,max=1.6;
  const svg=svgFor(document.getElementById("chartBlind"),W,H,"AlphaFold accuracy for PETase and 11 relatives; core RMSD under 0.45 angstroms for all");
  const x=v=>L+(v/max)*(W-L-Rr);
  [0,0.4,0.8,1.2,1.6].forEach(v=>{el("line",{x1:x(v),x2:x(v),y1:T,y2:B,stroke:css("--rule"),"stroke-width":1},svg);el("text",{x:x(v),y:B+16,"text-anchor":"middle"},svg).textContent=v.toFixed(1)});
  el("text",{x:(L+W-Rr)/2,y:B+32,"text-anchor":"middle"},svg).textContent="AlphaFold vs lab RMSD (Å)";
  let cutDrawn=false;
  rows.forEach((r,i)=>{
    const cy=T+i*rowH+rowH/2+4;
    if(r.blind&&!cutDrawn){cutDrawn=true;el("line",{x1:8,x2:W-Rr,y1:cy-rowH/2,y2:cy-rowH/2,stroke:css("--accent"),"stroke-width":1.5,"stroke-dasharray":"5 4"},svg);el("text",{x:W-Rr,y:cy-rowH/2-5,"text-anchor":"end",fill:css("--accent"),style:"fill:var(--accent)"},svg).textContent="solved after AlphaFold's training ↓"}
    el("text",{x:L-12,y:cy+4,"text-anchor":"end",class:"lbl"},svg).textContent=`${r.label} · ${r.year}`;
    if(r.strict)el("text",{x:x(Math.max(r.full,r.core))+12,y:cy+4,style:"fill:var(--accent);font-weight:600"},svg).textContent="strictly blind";
    el("line",{x1:x(r.core),x2:x(r.full),y1:cy,y2:cy,stroke:css("--s1"),"stroke-width":2,"stroke-opacity":.4},svg);
    const extra=(r.off?`<br>${r.off} residues >2 Å off, avg confidence ${r.offpl}`:"")+(r.mdate?`<br>model dated ${r.mdate}${r.strict?" (before its structure)":""}<br>core keeps ${r.pct}% of residues`:"");
    dot(svg,x(r.full),cy,5,css("--s1"),css("--bg"),`${r.label}<br>all residues ${f2(r.full)} Å<br>core ${f2(r.core)} Å${extra}`,true);
    dot(svg,x(r.core),cy,5,css("--s1"),css("--bg"),`${r.label}<br>core ${f2(r.core)} Å`);
  });
}

