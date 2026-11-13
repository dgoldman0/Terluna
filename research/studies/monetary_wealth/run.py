"""Reproduce annual monetary, wealth and ownership scenarios."""
from __future__ import annotations
from copy import deepcopy
import gzip,hashlib,json
from pathlib import Path
from . import model
HERE=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load():
    c=json.loads((HERE/'wealth_inputs.json').read_text())
    products={}
    for key in ('asset_product','source_access'):
        s=c[key];p=HERE/s['path']
        if sha(p)!=s['sha256']:raise ValueError('Changed '+key)
        x=json.loads(gzip.decompress(p.read_bytes()))
        if x['schema']!=s['schema']:raise ValueError('Unexpected '+key+' schema')
        products[key]=x
    p=HERE/'inputs/finance_flows.json.gz';f=json.loads(gzip.decompress(p.read_bytes()))
    if f['producer']['sha256']!=sha(HERE/'accounts.py'):raise ValueError('Stale finance flow producer')
    for name,expected in f['dependencies'].items():
        if sha(HERE/name)!=expected:raise ValueError('Changed finance dependency')
    source_finance=json.loads((HERE/'../enclosed_returns/finance_inputs.json').read_text())
    if c['owner_custody_fraction']!=source_finance['monetary']['incremental_distribution_spending_fraction']:
        raise ValueError('Custody allocation differs from the source finance policy')
    return c,f,products

def product():
    c,f,p=load();assets={a['name']:a for a in p['asset_product']['cases']};results=[];sens=[]
    for case in f['cases']:
        for ownership in c['ownership_scenarios']:
            results.append(model.path(case,assets[case['name']],c,ownership))
        if case['name']==c['reference_case']:
            for field,values in [('household_net_saving_fraction',c['household_saving_sensitivities']),('external_investor_U_retention_fraction',c['external_investor_U_retention_sensitivities']),('purchase_discount_rate',c['purchase_rate_sensitivities'])]:
                for value in values:
                    settings=deepcopy(c);settings[field]=value
                    for ownership in ('external_equity','community_purchase'):
                        result=model.path(case,assets[case['name']],settings,ownership)
                        sens.append({'parameter':field,'value':value,'ownership':ownership,'checkpoints':result['checkpoints']})
            result=model.path(case,assets[case['name']],c,'community_purchase',purchase_basis='own_funds')
            sens.append({'parameter':'negotiated_purchase_basis','value':'own_funds','ownership':'community_purchase','checkpoints':result['checkpoints']})
    return {'schema':'terluna.humanities.monetary-wealth/1','producer':{'path':'research/studies/monetary_wealth/run.py','sha256':sha(__file__)},'dependencies':{name:sha(HERE/name) for name in ('wealth_inputs.json','model.py','ledger.py','accounts.py','inputs/finance_flows.json.gz','inputs/asset_stock.json.gz','inputs/source_access.json.gz','sources.json','../complex_returns/consolidated.py')},'boundary':'Whole planned off-Earth enclosed zone, excluding Earth; no Open Moon, array or outer-system expansion. Outside investors/customers remain explicit counterparties.','evidence':'Conditional sector balance sheets and ownership/saving/monetary-holding scenarios on unchanged physical and financial paths; no calibrated money demand, nominal U exchange rate, legal ownership allocation or market valuation.','units':'constant USD2041 reference purchasing-power values under the U accounting convention','cases':results,'sensitivities':sens,'unvalued_assets':p['asset_product']['unvalued_assets']+['Complete atmospheric/water, consumable, cargo and product inventories are not exported; omitted values are unknown, not demonstrated zero.'],'rules':{k:v for k,v in c.items() if k.endswith('_rule') or k in ('currency_value_convention','household_closure','reserve_ownership')}}

def main():
    x=product();path=HERE/'results/wealth.json.gz';path.write_bytes(gzip.compress((json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0));print(path,len(x['cases']))

if __name__=='__main__':main()
