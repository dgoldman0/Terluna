"""Build the self-contained A2 accuracy lab from verified local evidence."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent;GEN=ROOT/'data/a2/generated'
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def build(output):
 val=load(ROOT/'A2_VALIDATION.json')
 if val['status']!='reference-checkpoint-passed':raise ValueError('A2 validation has not passed')
 sky=load(GEN/'moon-fair.multiple-sky.json');fine=load(GEN/'moon-fair-molecular-fine-full-chunked.json');cloud=load(GEN/'cloud-convergence.json');caps={i:load(GEN/f'capture-{i}.json') for i in ['grid-coarse','grid-a1','grid-fine','crop-075','crop-150','crop-200']}
 if sha(GEN/'moon-fair.multiple-sky.json')!=val['molecular']['atlas_sha256']:raise ValueError('Stale sky evidence')
 if sha(GEN/'moon-fair-molecular-fine-full-chunked.json')!=val['molecular']['fine_chunked_sha256']:raise ValueError('Stale fine evidence')
 data={'sky':sky,'fine':fine,'cloud':cloud,'captures':caps,'validation':val,'validationSha':sha(ROOT/'A2_VALIDATION.json')}
 payload=json.dumps(data,separators=(',',':')).replace('</','<\\/');html=(ROOT/'index.a2-accuracy-lab.html').read_text().replace('__A2_PAYLOAD__',payload).replace('__A2_SCRIPT__',(ROOT/'src/a2-accuracy-lab.js').read_text());output.write_text(html);print(json.dumps({'output':str(output),'bytes':output.stat().st_size,'sha256':sha(output)}))
if __name__=='__main__':build(ROOT/'Open_Moon_A2_Accuracy_Lab.html')
