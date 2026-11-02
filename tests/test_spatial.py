from tdln_pd1.spatial import assign_zone, zone_table


def test_zone_map_contains_lymph_node_microanatomy():
    assert {"subcapsular_sinus", "paracortex", "follicle", "medulla"}.issubset(set(zone_table().zone))
    assert assign_zone("CD8 T cell") == "paracortex"
    assert assign_zone("B cell") == "follicle"
