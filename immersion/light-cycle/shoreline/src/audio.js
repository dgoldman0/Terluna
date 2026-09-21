/* Procedural spatial audio. Acoustic textures are illustrative, not a sound-speed solver. */
(function(root){'use strict';const OM=root.OM,C=root.OpenMoonCore;
class Soundscape{
 constructor(){this.context=null;this.active=false;this.nodes=[];this.stepDistance=0;this.lastStep=0;}
 async start(){if(this.context){await this.context.resume();this.active=!this.active;return this.active;}const AC=root.AudioContext||root.webkitAudioContext;if(!AC)throw Error('Web Audio is unavailable.');const ctx=this.context=new AC();this.master=ctx.createGain();this.master.gain.value=.55;this.master.connect(ctx.destination);this.active=true;
  const rng=C.rng(19234),buffer=ctx.createBuffer(1,ctx.sampleRate*8,ctx.sampleRate),data=buffer.getChannelData(0);let pink=0;for(let i=0;i<data.length;i++){pink=.98*pink+.02*(rng()*2-1);data[i]=pink*3;}
  this.noise=buffer;
  const add=(name,freq,q,pos,amp)=>{const s=ctx.createBufferSource();s.buffer=buffer;s.loop=true;const f=ctx.createBiquadFilter();f.type='bandpass';f.frequency.value=freq;f.Q.value=q;const g=ctx.createGain();g.gain.value=amp;const p=ctx.createPanner();p.panningModel='HRTF';p.distanceModel='inverse';p.refDistance=12;p.maxDistance=250;p.rolloffFactor=1;p.positionX.value=pos[0];p.positionY.value=pos[1];p.positionZ.value=pos[2];s.connect(f).connect(g).connect(p).connect(this.master);s.start(0,(this.nodes.length*1.37)%7);const o={s,f,g,p};this.nodes.push(o);this[name]=o;};
  add('shore',480,.35,[0,1,-18],.4);add('leaves',1800,.3,[-18,8,40],.08);add('rain',4600,.22,[0,10,0],0);add('roof',850,.4,[20,8,33],0);add('drips',2400,1.5,[20,4,28],0);await ctx.resume();return true;
 }
 update(camera,w,ledger,time,velocity,muted=false){if(!this.context)return;const ctx=this.context,t=ctx.currentTime,roof=C.roofMask(camera.position.x,camera.position.z);this.master.gain.setTargetAtTime(this.active&&!muted?.55:0,t,.2);
  const n=this.shore;const surge=.35+.65*(.5+.5*Math.sin(time*.22+Math.sin(time*.041)));n.g.gain.setTargetAtTime((.10+w.wind*.023)*surge,t,.3);
  this.leaves.g.gain.setTargetAtTime((.015+w.wind*.018)*(1-.3*roof),t,.2);this.leaves.f.frequency.setTargetAtTime(900+w.wind*160,t,.3);
  this.rain.p.positionX.value=camera.position.x;this.rain.p.positionY.value=camera.position.y+6;this.rain.p.positionZ.value=camera.position.z;this.rain.g.gain.setTargetAtTime(w.rain*.022*(roof?.17:1),t,.08);this.roof.g.gain.setTargetAtTime(w.rain*.034,t,.1);this.drips.g.gain.setTargetAtTime(ledger.leaf*.12*(.55+.45*Math.sin(time*2.3)**6),t,.04);
  const listener=ctx.listener,dir=camera.getWorldDirection(new root.THREE.Vector3());for(const[k,v]of Object.entries({positionX:camera.position.x,positionY:camera.position.y,positionZ:camera.position.z,forwardX:dir.x,forwardY:dir.y,forwardZ:dir.z,upX:0,upY:1,upZ:0})){if(listener[k])listener[k].setTargetAtTime(v,t,.03);}
  if(velocity>.12&&t-this.lastStep>.53&&this.active&&!muted){this.lastStep=t;this.footstep(roof,ledger.exposed);}
 }
 footstep(wood,wet){const ctx=this.context,src=ctx.createBufferSource();src.buffer=this.noise;const f=ctx.createBiquadFilter();f.type='lowpass';f.frequency.value=wood?850:1700;const g=ctx.createGain(),t=ctx.currentTime;g.gain.setValueAtTime(.0001,t);g.gain.exponentialRampToValueAtTime(wood?.55:.42,t+.008);g.gain.exponentialRampToValueAtTime(.0001,t+.14+Math.min(wet,.8)*.07);src.connect(f).connect(g).connect(this.master);src.start(t,Math.random()*5,.23);src.onended=()=>{src.disconnect();f.disconnect();g.disconnect();};}
}
OM.Soundscape=Soundscape;
})(globalThis);
