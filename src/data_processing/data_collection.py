#!/usr/bin/env python3
"""
Step 02: Enhanced Web Data Collection from ClinicalTrials.gov API (v2)

- Queries studies (structured metadata + interventions)
- Extracts intervention names/types using armsInterventionsModule
- Outputs raw trial data to CSV including drug/device names in structured format
"""

import requests
import pandas as pd
from pathlib import Path
import time
import json
from tqdm import tqdm


class ClinicalTrialsWebCollector:
    def __init__(self, output_dir="data/raw"):
        self.base_url = "https://clinicaltrials.gov/api/v2/studies"
        self.rate_limit = 1.2  # seconds between requests (politeness delay)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def search_trials(self, query_term="", max_studies=500):
        all_studies, page_token, page_size = [], None, 1000
        pbar = tqdm(total=max_studies, desc=f"🔍 Fetching '{query_term}' studies")

        while len(all_studies) < max_studies:
            try:
                params = {
                    "format": "json",
                    "pageSize": min(page_size, max_studies - len(all_studies))
                }
                if query_term:
                    params["query.cond"] = query_term
                if page_token:
                    params["pageToken"] = page_token

                r = requests.get(self.base_url, params=params, timeout=30)
                r.raise_for_status()
                data = r.json()

                studies = data.get("studies", [])
                if not studies:
                    break

                processed = [self._extract_fields(study) for study in studies]
                processed = [p for p in processed if p]  # remove None

                all_studies.extend(processed)
                pbar.update(len(processed))

                page_token = data.get("nextPageToken")
                if not page_token:
                    break

                time.sleep(self.rate_limit)

            except Exception as e:
                print(f"⚠️ Error during API call: {e}")
                break

        pbar.close()

        df = pd.DataFrame(all_studies[:max_studies])
        if not df.empty:
            timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
            out_file = self.output_dir / f"raw_clinical_trials_{query_term or 'all'}_{timestamp}.csv"
            df.to_csv(out_file, index=False)
            print(f"✅ Saved {len(df)} trials to {out_file}")
        else:
            print("⚠️ No trials retrieved from ClinicalTrials.gov")

                df = pd.DataFrame(all_studies[:max_studies])
        if not df.empty:
            timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
            out_file = self.output_dir / f"raw_clinical_trials_{query_term or 'all'}_{timestamp}.csv"
            df.to_csv(out_file, index=False)
            print(f"✅ Saved {len(df)} trials to {out_file}")

            # ✅ Extract and save unique drug/intervention names
            if "interventions_names" in df.columns:
                interventions_col = df["interventions_names"].fillna("")
                intervention_names = set()
                for row in interventions_col:
                    for name in row.split(";"):
                        clean_name = name.strip().lower()
                        if clean_name:
                            intervention_names.add(clean_name)

                drugs_df = pd.DataFrame({"drug_name": sorted(intervention_names)})
                drugs_df.to_csv("data/processed/unique_drugs.csv", index=False)
                print(f"✅ Extracted and saved {len(drugs_df)} unique drugs to: data/processed/unique_drugs.csv")
            else:
                print("⚠️ Column 'interventions_names' not found in trials data.")

        else:
            print("⚠️ No trials retrieved from ClinicalTrials.gov")
        return df

    def _extract_fields(self, study):
        try:
            protocol = study.get('protocolSection', {})
            identification = protocol.get('identificationModule', {})
            description = protocol.get('descriptionModule', {})
            status = protocol.get('statusModule', {})
            design = protocol.get('designModule', {})
            design_info = design.get('designInfo', {})
            masking_info = design_info.get('maskingInfo', {})
            enrollment_info = design.get('enrollmentInfo', {})
            arms_module = protocol.get('armsInterventionsModule', {})

            interventions_list = []
            for intervention in arms_module.get('interventions', []):
                interventions_list.append({
                    "name": intervention.get("name", ""),
                    "type": intervention.get("type", ""),
                    "description": intervention.get("description", ""),
                    "mesh_terms": intervention.get("meshTerms", []),
                    "other_ids": intervention.get("otherIds", []),
                })

            interventions_json = json.dumps(interventions_list, ensure_ascii=False)

            return {
                "nct_id": identification.get("nctId", ""),
                "brief_title": identification.get("briefTitle", ""),
                "official_title": identification.get("officialTitle", ""),
                "brief_summary": description.get("briefSummary", ""),
                "detailed_description": description.get("detailedDescription", ""),
                "overall_status": status.get("overallStatus", ""),
                "why_stopped": status.get("whyStopped", ""),
                "start_date": status.get("startDateStruct", {}).get("date", ""),
                "completion_date": status.get("completionDateStruct", {}).get("date", ""),
                "primary_completion_date": status.get("primaryCompletionDateStruct", {}).get("date", ""),
                "study_type": design.get("studyType", ""),
                "phases": '|'.join(design.get("phases", [])) if design.get("phases") else "",
                "allocation": design_info.get("allocation", ""),
                "intervention_model": design_info.get("interventionModel", ""),
                "masking": masking_info.get("masking", ""),
                "primary_purpose": design_info.get("primaryPurpose", ""),
                "enrollment_count": enrollment_info.get("count", 0),
                "enrollment_type": enrollment_info.get("type", ""),
                "interventions_json": interventions_json,
                "interventions_types": ";".join(i.get("type", "") for i in interventions_list if i.get("type", "")),
                "interventions_names": ";".join(i.get("name", "") for i in interventions_list if i.get("name", ""))
            }

        except Exception as e:
            print(f"❌ Extraction error: {e}")
            return None


if __name__ == "__main__":
    collector = ClinicalTrialsWebCollector(output_dir="data/raw")

    # Run multiple queries: can adjust or remove as needed
    collector.search_trials("cancer", max_studies=500)
    collector.search_trials("diabetes", max_studies=500)
    collector.search_trials("hypertension", max_studies=500)

    print("🚀 ClinicalTrials.gov data collection complete.")
