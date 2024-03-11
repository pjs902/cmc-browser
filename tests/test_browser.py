from cmcbrowser import CMCBrowser
from tqdm import tqdm


def test_all_snapshot_types():
    b = CMCBrowser(ss_dir="/home/peter/research/CMC-grid/grid/")

    print("Found the following models:")
    print(b.models_list)

    for model in tqdm(b.models_list):
        print(f"Testing regular snapshots for model: {model}")

        snaps_regular = b.model_snapshots[model]["regular"]
        for ss in snaps_regular:
            if "dat.gz" in ss:
                b.load_snapshot(
                    model_name=model,
                    ss_name=ss,
                    distance=5,
                    mode="dat.gz",
                    strict=False,
                )

                snap = b.loaded_snapshots[f"{model}/{ss}"]
                print("Successfully loaded snapshot")
                break

            else:
                b.load_snapshot(
                    model_name=model, ss_name="king.snapshots.h5", distance=5, mode="h5", strict=False, h5_key=ss
                )
                snap = b.loaded_snapshots[f"{model}/{ss}"]
                print("Successfully loaded snapshot")
                break

        print("Testing window snapshots for model: {model}")
        snaps_window = b.model_snapshots[model]["window"]

        if snaps_window is None:
            print(f"No window snapshots for {model}")
            continue

        for ss in snaps_window:
            b.load_snapshot(
                model_name=model, ss_name="king.window.snapshots.h5", distance=5, mode="h5", strict=False, h5_key=ss
            )
            snap = b.loaded_snapshots[f"{model}/{ss}"]
            print("Successfully loaded snapshot")
            break


if __name__ == "__main__":
    test_all_snapshot_types()
