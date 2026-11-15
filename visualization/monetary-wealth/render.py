"""Render the conditional money and wealth accounts without changing them."""
from pathlib import Path
import gzip, hashlib, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'research/studies/monetary_wealth/results/wealth.json.gz'
OUT=ROOT/'research/runs/monetary_wealth'

def main():
    product=json.loads(gzip.decompress(SOURCE.read_bytes()))
    cases={c['ownership']:c for c in product['cases'] if c['name']=='propellant_logistics:breakthrough:wider'}
    baseline=cases['external_equity']['annual'];years=[r['year'] for r in baseline]
    series={
        'gross_operating_assets':[r['physical_assets']['gross_replacement_value_USD2041']/1e12 for r in baseline],
        'net_assets_and_construction':[r['book_balance']['consolidated_off_earth']['real_assets']/1e12 for r in baseline],
        'enterprise_own_funds':[r['book_balance']['own_funds_equity']/1e12 for r in baseline]}
    fig,axes=plt.subplots(2,2,figsize=(13.5,9.5))
    for key,label,color in [('gross_operating_assets','Gross operating replacement cost','#557a95'),('net_assets_and_construction','Remaining service + construction','#195b43'),('enterprise_own_funds','Enterprise own funds','#9b6728')]:
        axes[0,0].plot(years,series[key],label=label,color=color,lw=2)
    for ownership,label,color in [('external_equity','Continued outside ownership','#596579'),('community_purchase','Paid community purchases, 3% quote','#dc8532'),('negotiated_vesting','Negotiated transfer to 75%','#237e67')]:
        rows=cases[ownership]['annual']
        for field,axis,values in [('transaction_money',axes[0,1],[r['money_stock_flow']['outstanding_domestic_issuer_nonbank_money']/1e12 for r in rows]),('local_net_wealth',axes[1,0],[r['book_balance']['consolidated_off_earth']['net_worth']/1e12 for r in rows]),('outside_cash_receipts',axes[1,1],[(r['cumulative']['outside_payouts']+r['cumulative']['outside_buyout_receipts'])/1e12 for r in rows])]:
            series[ownership+':'+field]=values
            axis.plot(years,values,label=label,color=color,lw=2)
    titles=['Productive assets and enterprise claims','Outstanding U transaction-money value','Resident and community net wealth — logarithmic scale','Cumulative outside distributions and sale receipts']
    for ax,title in zip(axes.flat,titles):
        ax.set_title(title,loc='left',fontsize=11,fontweight='bold',pad=10)
        ax.set_xlabel('Year from initial funding');ax.set_ylabel('Trillion USD2041 equivalent')
        ax.set_xlim(0,220);ax.set_xticks([0,32,100,200,220]);ax.grid(alpha=.18)
        for year in (100,200):ax.axvline(year,color='#666666',alpha=.25,lw=.8)
        ax.legend(fontsize=8,frameon=False,loc='upper left')
    axes[1,0].set_yscale('log');axes[1,0].set_ylim(1e-6,1000)
    fig.suptitle('U monetary balances and ownership through year 220\nPropellant and worldwide services · breakthrough demand · 1% real growth to year 200',fontsize=15,fontweight='bold')
    fig.subplots_adjust(left=.075,right=.98,bottom=.13,top=.86,hspace=.34,wspace=.22)
    fig.text(.02,.025,'Conditional accounts: 15% payroll saving; zero real yield; community retains all payouts. Transfers require investor agreement.\nMoney is a reference-valued stock, not nominal U units. Physical supply and deployment remain unqualified.',fontsize=10)
    OUT.mkdir(parents=True,exist_ok=True)
    for ext in ('png','svg'):fig.savefig(OUT/('wealth.'+ext),dpi=160)
    manifest={'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'units':'trillion USD2041 reference purchasing-power equivalent','years':years,'series':series}
    (OUT/'wealth.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(OUT/'wealth.png')
if __name__=='__main__':main()
