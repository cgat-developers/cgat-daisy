import pandas
import sqlalchemy
import matplotlib.pyplot as plt

database = sqlalchemy.create_engine("sqlite:///./csvdb")

df = pandas.read_sql(
    "SELECT i.tool_name, m.key, m.value "
    "FROM bcftools_stats_summary_numbers AS m, instance AS i "
    "WHERE i.id = m.instance_id AND key='number_of_SNPs' ", database).set_index("tool_name")

df.plot.bar()
plt.tight_layout()
plt.savefig("number_variants.png")
