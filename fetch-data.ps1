# Downloads every input used by compare.ps1, ensemble.ps1 and family.ps1 into data\.
# Safe to re-run: existing files are kept.
$Root = $PSScriptRoot; $D = Join-Path $Root "data"
foreach ($sub in "", "wt", "family") { New-Item -ItemType Directory -Force (Join-Path $D $sub) | Out-Null }
function Get-File($url, $path) { if (-not (Test-Path $path)) { Invoke-WebRequest $url -OutFile $path -UseBasicParsing } }

# PETase sequence, AlphaFold model, and the single-comparison structures
Get-File "https://rest.uniprot.org/uniprotkb/A0A0K8P6T7.fasta" (Join-Path $D "PETase_A0A0K8P6T7.fasta")
$af = Invoke-RestMethod "https://alphafold.ebi.ac.uk/api/prediction/A0A0K8P6T7"
Get-File $af[0].pdbUrl (Join-Path $D "AF-A0A0K8P6T7_predicted.pdb")
foreach ($id in "5XJH", "6EQE") { Get-File "https://files.rcsb.org/download/$id.pdb" (Join-Path $D "${id}_experimental.pdb") }
Get-File "https://files.rcsb.org/download/7SH6.pdb" (Join-Path $D "7SH6_fastpetase.pdb")   # FAST-PETase, for part2.ps1

function Search-Entities($attribute, $value) {
  $q = @{ query = @{ type = "terminal"; service = "text"; parameters = @{ attribute = $attribute; operator = "exact_match"; value = $value } }
          return_type = "polymer_entity"; request_options = @{ return_all_hits = $true } } | ConvertTo-Json -Depth 6 -Compress
  @((Invoke-RestMethod -Method Post -Uri "https://search.rcsb.org/rcsbsearch/v2/query" -Body $q -ContentType "application/json").result_set | ForEach-Object identifier)
}
function Get-Metadata($ids) {
  $out = @()
  for ($i = 0; $i -lt $ids.Count; $i += 100) {
    $chunk = $ids[$i..([math]::Min($i + 99, $ids.Count - 1))] | ForEach-Object { '"' + $_ + '"' }
    $gql = 'query { polymer_entities(entity_ids:[' + ($chunk -join ',') + ']) { rcsb_id entity_poly { rcsb_mutation_count } rcsb_polymer_entity_container_identifiers { uniprot_ids auth_asym_ids } rcsb_entity_source_organism { scientific_name } rcsb_polymer_entity { pdbx_description } entry { rcsb_accession_info { initial_release_date } rcsb_entry_info { resolution_combined experimental_method nonpolymer_bound_components } struct { title } } } }'
    $r = Invoke-RestMethod -Method Post -Uri "https://data.rcsb.org/graphql" -Body (@{ query = $gql } | ConvertTo-Json -Compress) -ContentType "application/json"
    $out += foreach ($e in $r.data.polymer_entities) { [pscustomobject]@{
      entity = $e.rcsb_id; id = $e.rcsb_id.Split('_')[0]; uniprot = ($e.rcsb_polymer_entity_container_identifiers.uniprot_ids -join ',')
      muts = $e.entity_poly.rcsb_mutation_count; org = ($e.rcsb_entity_source_organism.scientific_name | Select-Object -First 1)
      desc = $e.rcsb_polymer_entity.pdbx_description; released = $e.entry.rcsb_accession_info.initial_release_date.Substring(0, 10)
      res = ($e.entry.rcsb_entry_info.resolution_combined | Select-Object -First 1); method = $e.entry.rcsb_entry_info.experimental_method
      chains = ($e.rcsb_polymer_entity_container_identifiers.auth_asym_ids -join ''); ligands = ($e.entry.rcsb_entry_info.nonpolymer_bound_components -join ',')
      title = $e.entry.struct.title } }
  }
  $out
}

# Every PETase structure. Candidates are those the database marks unmutated; ensemble.ps1 then
# checks the actual sequences and drops the ones that aren't wild-type.
$pet = Get-Metadata (Search-Entities "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession" "A0A0K8P6T7")
$pet | Export-Csv (Join-Path $D "petase_structures.csv") -NoTypeInformation -Encoding UTF8
$pet | Where-Object { $_.muts -eq 0 } | ForEach-Object { Get-File "https://files.rcsb.org/download/$($_.id).pdb" (Join-Path $D "wt\$($_.id).pdb") }

# PET hydrolase family (EC 3.1.1.101); family.ps1 downloads each relative's structure and model
$fam = Get-Metadata (Search-Entities "rcsb_polymer_entity.rcsb_ec_lineage.id" "3.1.1.101")
$fam | Export-Csv (Join-Path $D "pet_hydrolase_family.csv") -NoTypeInformation -Encoding UTF8
"Done: $($pet.Count) PETase entities, $($fam.Count) family entities."
