from tqdm import tqdm


def extract_ca_atoms(input_file, output_file, res_file, step, model_start=1, model_end=50000, special_res=[]):
    # 计算总行数
    with open(input_file, "r") as f:
        total_lines = sum(1 for _ in f)

    with open(input_file, "r") as f_in, open(output_file, "w") as f_out, open(res_file, "w") as f_res:
        pass_list = ["REMARK", "TITLE", "CRYST1", "TER"]
        atom_count = 1
        res_ls = []
        progress_bar = tqdm(total=total_lines, desc="Processing", unit="line")

        for line in f_in:
            progress_bar.update(1)

            if any(line.startswith(keyword) for keyword in pass_list):
                continue
            elif line.startswith("ATOM"):
                if atom_count % step == 0:
                    res_before = int(line[4:11])
                    line = line[:6] + f"{int(atom_count / step):>5d}" + line[
                                                                        11:21] + f"{int(atom_count / step):>5d}" + line[
                                                                                                                   26:]
                    res_alter = int(line[4:11])
                    res_line = str(res_before) + "  ----->  " + str(res_alter) + "\n"
                    if res_line not in res_ls:
                        res_ls.append(res_line)
                    f_out.write(line)
                elif int(line[23:26]) in special_res:
                    atom_count += step
                    res_before = int(line[4:11])
                    line = line[:6] + f"{int(atom_count / step):>5d}" + line[
                                                                        11:21] + f"{int(atom_count / step):>5d}" + line[
                                                                                                                   26:]
                    res_alter = int(line[4:11])
                    res_line = str(res_before) + "  ----->  " + str(res_alter) + "\n"
                    if res_line not in res_ls:
                        res_ls.append(res_line)
                    f_out.write(line)
                atom_count += 1
            elif line.startswith("MODEL"):
                if int(line[8:]) < model_start:
                    continue
                elif int(line[8:]) > model_end:
                    break
                else:
                    atom_count = 1
                    f_out.write(line)
            else:
                f_out.write(line)

        for res in res_ls:
            f_res.write(res)
        progress_bar.close()



extract_ca_atoms(f"/home/ps/NRI-sir/ca_1.pdb",
                 f"/home/ps/NRI-sir/data/pdb/ca_1.pdb",
                 r"/home/ps/NRI-sir/res_index.txt",
                 2, 1, 5000)
                     
