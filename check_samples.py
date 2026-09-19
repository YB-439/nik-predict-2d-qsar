import json

data = json.load(open("coralsea_full_predictions.json"))
print("First 10 dataset compounds:")
for k, v in list(data.items())[:10]:
    print(f"{v['code']}: exp={v['exp_pIC50']}, pred={v['consensus']}, smiles={k}")
