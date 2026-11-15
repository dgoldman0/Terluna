"""Export unchanged dated finance flows for a sector/ownership consumer."""
from __future__ import annotations
import gzip,hashlib,json
from pathlib import Path
from ..enclosed_returns import finance
from ..complex_returns import credit
HERE=Path(__file__).resolve().parent

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def product():
    cfg=json.loads((HERE/'wealth_inputs.json').read_text())
    for name,expected in cfg['finance_dependencies'].items():
        if digest(HERE/name)!=expected:raise ValueError('Changed finance input '+name)
    fc,orders,products=finance.inputs();requests={c['name']:c for c in orders['cases']}
    checkpoints=json.loads((HERE/'../enclosed_returns/results/returns.json').read_text())['checkpoints'];result=[]
    for physical in products['physical']['cases']:
        if physical.get('physical_scenario','base')!='base' or physical['market']!=cfg['market'] or physical['market_growth_name'] not in cfg['growth_cases']:continue
        name=physical['name'];p=finance.plan(requests[name],physical,fc);control=finance.finance(p,fc,'equity');f=finance.finance(p,fc,'U_credit',control)
        for horizon in (100,200,220):
            computed=finance.metrics(p,f,fc,horizon);old=next(c for c in checkpoints if c['name']==name and c['horizon']==horizon and c['mode']=='U_credit')
            for key in ('owner_payouts','terminal_cash','terminal_reserves','terminal_debt','terminal_custody','additional_equity_required'):
                if abs(computed[key]-old[key])>max(1.,abs(old[key])*1e-12):raise ValueError('Baseline changed '+name+' '+key)
        rows=[];loans=[];start=p['annual'][0]['year']
        for y in range(start):
            funding=p['seed_hardware_and_carry']/start+(f['initial_reserve'] if y==0 else 0)+(p['initial_working'] if y==start-1 else 0)
            rows.append({'year':y,'population_convention':p['seed_population']/2,'gross':0.,'owner_payout':0.,'dividend':0.,'owner_custody':0.,'resident_payroll':p['seed_resident_payroll']/start,'cash':p['initial_working'] if y==start-1 else 0.,'reserves':f['initial_reserve'],'outside_equity_contribution':funding,'contribution':0.,'debt':0.,'local_U_debt':0.,'supplier_U_debt':0.,'external_debt':0.,'local_U_draw':0.,'supplier_U_draw':0.,'external_draw':0.,'local_U_principal':0.,'supplier_U_principal':0.,'external_principal':0.,'local_U_interest':0.,'supplier_interest':0.,'external_interest':0.,'resource_provision':physical['seed']['total_resource_provision_USD2041']/start,'external_operating':None,'capital':None,'service_residual_capital':None,'service_residual_renewal':None})
        for path,flow in zip(p['annual'],f['annual']):
            y=flow['year'];due=credit.repay(loans,y)
            for field,kind,rate,share in [('local_U_draw','local',fc['local_U_rate'],0),('supplier_U_draw','supplier',fc['external_rate'],1),('external_draw','external',fc['external_rate'],0)]:
                draw=flow[field]
                if draw>1e-6:loans.append({'year':y,'original':draw,'remaining':draw,'rate':rate,'term':fc['loan_term_years'],'grace':fc['underwriting']['principal_grace_years'],'kind':'external' if kind=='external' else 'U','supplier_share':share,'sector_kind':kind})
            balances={k:sum(l['remaining'] for l in loans if l['sector_kind']==k) for k in ('local','supplier','external')}
            if abs(sum(balances.values())-flow['debt'])>1:raise ValueError('Debt cohort mismatch')
            rows.append({**path,**flow,'outside_equity_contribution':flow['contribution'],'local_U_debt':balances['local'],'supplier_U_debt':balances['supplier'],'external_debt':balances['external'],'local_U_principal':due['U_principal']-due['supplier_principal'],'supplier_U_principal':due['supplier_principal'],'external_principal':due['principal']-due['U_principal']})
        result.append({'name':name,'portfolio':p['portfolio'],'growth':p['growth'],'market':p['market'],'seed':p['seed'],'initial_reserve':f['initial_reserve'],'seed_cost_scope':p['seed_cost_scope'],'unpriced_deployment_requirements':p['unpriced_deployment_requirements'],'annual':rows,'baseline_checkpoints':[c for c in checkpoints if c['name']==name and c['mode']=='U_credit']})
    if len(result)!=18:raise ValueError('Expected 18 portfolio/growth paths')
    return {'schema':'terluna.humanities.wealth-finance-flows/1','producer':{'path':'research/studies/monetary_wealth/accounts.py','sha256':digest(__file__)},'dependencies':{**cfg['finance_dependencies'],'../enclosed_returns/results/orders.json':digest(HERE/'../enclosed_returns/results/orders.json'),'../enclosed_returns/inputs/industry.json.gz':digest(HERE/'../enclosed_returns/inputs/industry.json.gz')},'cases':result,'evidence':'Unchanged finance model rerun; all 54 saved checkpoint cash/reserve/debt/funding/payout comparisons reproduce. Complete annual 0–220 ledger with seed cash timing inherited from finance.metrics.'}

def main():
    x=product();p=HERE/'inputs/finance_flows.json.gz';p.write_bytes(gzip.compress((json.dumps(x,separators=(',',':'),sort_keys=True)+'\n').encode(),mtime=0));print(p,len(x['cases']))

if __name__=='__main__':main()
