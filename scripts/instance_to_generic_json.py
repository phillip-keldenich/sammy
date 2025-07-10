import itertools, json, re, sys
from samplns.instances import (
    parse as parse_xml_instance,
    FeatureNode,
    CompositeFeature,
    VAR,
    OR,
)
from samplns.instances.feature import (
    AltFeature,
    AndFeature,
    ConcreteFeature,
    FeatureLiteral,
    OrFeature,
)
from samplns.preprocessor import Preprocessor, IndexInstance


def _to_sat_literal_from_value(var_name: int, value: bool) -> int:
    v = var_name + 1
    if not value:
        v = -v
    return v


def _to_sat_literal(literal: FeatureLiteral | VAR) -> int:
    return _to_sat_literal_from_value(literal.var_name, not literal.negated)


class ClauseListModel:
    def __init__(self):
        self.clauses = []

    def add_clause(self, clause: list[int]):
        self.clauses.append(clause.copy())


class ClauseModelCreator:
    def __init__(self, model_object):
        self.model_object = model_object

    def __add_and_feature_constraint(self, node: AndFeature):
        par = node.feature_literal
        for e in node.elements:
            if not e.mandatory:
                cld = e.feature_literal
                self.model_object.add_clause(
                    [-_to_sat_literal(cld), _to_sat_literal(par)]
                )
            else:
                assert e.feature_literal.var_name == par.var_name
            self.__add_structure_constraints(e)

    def __add_or_feature_constraint(self, node: OrFeature):
        par = node.feature_literal
        lits = [_to_sat_literal(e.feature_literal) for e in node.elements]
        lits.append(-_to_sat_literal(par))
        self.model_object.add_clause(lits)
        for e in node.elements:
            cld = e.feature_literal
            self.model_object.add_clause([-_to_sat_literal(cld), _to_sat_literal(par)])
            self.__add_structure_constraints(e)

    def __add_alt_feature_constraint(self, node: AltFeature):
        par = node.feature_literal
        lits = [_to_sat_literal(e.feature_literal) for e in node.elements]
        for l1, l2 in itertools.combinations(lits, 2):
            self.model_object.add_clause([-l1, -l2])
        lits.append(-_to_sat_literal(par))
        self.model_object.add_clause(lits)
        for e in node.elements:
            cld = e.feature_literal
            self.model_object.add_clause([-_to_sat_literal(cld), _to_sat_literal(par)])
            self.__add_structure_constraints(e)

    def __add_structure_constraints(self, node):
        if isinstance(node, ConcreteFeature):
            return
        if isinstance(node, AndFeature):
            self.__add_and_feature_constraint(node)
        elif isinstance(node, OrFeature):
            self.__add_or_feature_constraint(node)
        elif isinstance(node, AltFeature):
            self.__add_alt_feature_constraint(node)
        else:
            raise ValueError("Unexpected ndoe in __add_structure_constraints!")

    def __add_rule_constraint(self, rule):
        if isinstance(rule, VAR):
            self.model_object.add_clause([_to_sat_literal(rule)])
        elif isinstance(rule, OR):
            self.model_object.add_clause([_to_sat_literal(e) for e in rule.elements])
        else:
            raise ValueError(f"Fomula not in CNF: {rule}")

    def create(self, instance: IndexInstance):
        self.model_object.add_clause(
            [_to_sat_literal(instance.structure.feature_literal)]
        )
        self.__add_structure_constraints(instance.structure)
        for rule in instance.rules:
            self.__add_rule_constraint(rule)


def create_clause_list(instance: IndexInstance):
    model = ClauseListModel()
    cmc = ClauseModelCreator(model)
    cmc.create(instance)
    return model.clauses


def _sort_structure(node: FeatureNode):
    if isinstance(node, CompositeFeature):
        for e in node.elements:
            _sort_structure(e)
        node.elements.sort(
            key=lambda x: (x.feature_literal.var_name, x.feature_literal.negated)
        )


class DIMACSParser:
    _p_line = re.compile(r"^p\s+cnf\s+([0-9]+)\s+([0-9]+)\s*$")
    _c_line = re.compile(r"^c\s+([0-9]+)(\$?)\s+([()'A-Za-z0-9_.,-]+)$")

    def __init__(self, file_path):
        with open(file_path, "r") as fstream:
            self.lines = fstream.readlines()
        self.clauses = []
        self.num_vars, self.num_clauses = None, None
        self.concrete_vars = []
        self.current_clause = []

    def _handle_comment(self, l):
        m = DIMACSParser._c_line.match(l)
        if not m:
            print("Warning: Unrecognized comment line", l)
        else:
            if len(m[2]) == 0:
                var_id = int(m[1])
                if var_id != len(self.concrete_vars) + 1:
                    raise ValueError("Concrete variables out of order!")
                self.concrete_vars.append(m[3])

    def _handle_formula_line(self, l):
        if self.num_vars is not None:
            raise ValueError("Multiple formula lines (p ...)")
        m = DIMACSParser._p_line.match(l)
        if not m:
            raise ValueError(f"Invalid formula line ({l})!")
        self.num_vars = int(m[1])
        self.num_clauses = int(m[2])

    def _handle_clause_line(self, l):
        split = l.split(" ")
        for iword in split:
            if len(iword) == 0:
                continue
            value = int(iword)
            if value == 0:
                self.clauses.append(self.current_clause)
                self.current_clause = []
            else:
                if abs(value) > self.num_vars:
                    raise ValueError(
                        f"Invalid literal (out of range): {value} (only {self.num_vars} variables)!"
                    )
                self.current_clause.append(value)

    def _handle_line(self, l):
        if l.startswith("c"):
            self._handle_comment(l)
        elif l.startswith("p"):
            self._handle_formula_line(l)
        else:
            if self.num_vars is None:
                raise ValueError("Expected formula line before first clause!")
            self._handle_clause_line(l)

    def parse(self):
        for l in self.lines:
            self._handle_line(l.strip())
        if len(self.current_clause) != 0:
            raise ValueError("Unfinished clause in file!")
        if len(self.clauses) != self.num_clauses:
            raise ValueError("Number of clauses does not match formula declaration!")


class GenericJSONInstance:
    def __init__(self, name):
        self.type: str = "software configuration model"
        self.meta: dict[str, object] = {}
        self.name: str = name
        self.labels: dict[str, int] = dict()
        self.num_variables: int = 0
        self.num_concrete_features: int = 0
        self.cnf_clauses: list[list[int]] = []
        self.structure: dict[str, object] | None = None

    def from_commented_dimacs(self, dimacs_file):
        parser = DIMACSParser(dimacs_file)
        parser.parse()
        self.num_variables = parser.num_vars
        self.num_concrete_features = len(parser.concrete_vars)
        self.cnf_clauses = parser.clauses
        self.structure = None
        for index, feature in enumerate(parser.concrete_vars):
            self.labels[feature] = index + 1

    def from_xml_model(self, xml_or_tgz):
        instance = parse_xml_instance(xml_or_tgz)
        instance.features.sort()
        _sort_structure(instance.structure)
        index_instance = Preprocessor(cnf=True).preprocess(instance)
        self.cnf_clauses = create_clause_list(index_instance)
        self.num_variables = index_instance.n_all
        self.num_concrete_features = index_instance.n_concrete
        self.structure = {
            "tree": instance.structure.to_json_data(),
            "rules": [rule.to_json_data() for rule in instance.rules],
        }
        for feature in instance.features:
            result = index_instance.to_mapped_universe({feature: True}, 
                                                       strict=False)
            assert len(result) == 1
            for s, v in result.items():
                l = abs(s) + 1
                self.labels[feature] = l if v else -l

    def __write(self, out_stream, **kwargs):
        if self.num_variables <= 0 or self.num_concrete_features <= 0:
            raise ValueError(f"Invalid instance: no variables or no concrete features!")
        json_data = {
            "type": self.type,
            "name": self.name,
            "meta": self.meta,
            "num_variables": self.num_variables,
            "num_concrete_features": self.num_concrete_features,
            "labels": self.labels,
            "cnf_clauses": self.cnf_clauses,
        }
        if self.structure:
            json_data["structure"] = self.structure
        json.dump(json_data, out_stream, **kwargs)

    def write(self, out_file, **kwargs):
        try:
            with open(out_file, "w") as f:
                self.__write(f, **kwargs)
        except TypeError:
            self.__write(out_file, **kwargs)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        input_file = sys.argv[1]
        if input_file.endswith(".xml"):
            output_file = input_file[:-4] + ".scm.json"
        elif input_file.endswith(".dimacs"):
            output_file = input_file[:-7] + ".scm.json"
        else:
            print("Expected either XML or DIMACS file!")
            exit(1)
    elif len(sys.argv) != 3:
        print("Expecting at least one and at most two parameters: <input file> [output file]")
        exit(1)
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
    if input_file.endswith("xml"):
        instance = GenericJSONInstance(input_file)
        instance.from_xml_model(input_file)
    else:
        instance = GenericJSONInstance(input_file)
        instance.from_commented_dimacs(input_file)
    instance.write(output_file, indent=2)
