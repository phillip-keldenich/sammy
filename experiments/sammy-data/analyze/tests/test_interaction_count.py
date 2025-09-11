from sample_analyzer.verify import count_interactions_py, count_interactions

sample = [
  {
    "Addressbook": True,
    "AutoResponder": True,
    "Base": True,
    "Decrypt": True,
    "Encrypt": True,
    "Forward": True,
    "Keys": True,
    "Sign": True,
    "Verify": True
  },
  {
    "Addressbook": True,
    "AutoResponder": True,
    "Base": True,
    "Decrypt": True,
    "Encrypt": True,
    "Forward": False,
    "Keys": True,
    "Sign": True,
    "Verify": True
  },
  {
    "Addressbook": False,
    "AutoResponder": True,
    "Base": True,
    "Decrypt": False,
    "Encrypt": False,
    "Forward": True,
    "Keys": False,
    "Sign": False,
    "Verify": False
  },
  {
    "Addressbook": True,
    "AutoResponder": False,
    "Base": True,
    "Decrypt": False,
    "Encrypt": False,
    "Forward": False,
    "Keys": False,
    "Sign": False,
    "Verify": False
  },
  {
    "Addressbook": False,
    "AutoResponder": False,
    "Base": True,
    "Decrypt": False,
    "Encrypt": False,
    "Forward": True,
    "Keys": True,
    "Sign": True,
    "Verify": True
  },
  {
    "Addressbook": False,
    "AutoResponder": False,
    "Base": True,
    "Decrypt": True,
    "Encrypt": True,
    "Forward": False,
    "Keys": True,
    "Sign": False,
    "Verify": False
  }
]


def test_interaction_count() -> int:
    concrete_features = [
        "Base",
        "Keys",
        "Encrypt",
        "AutoResponder",
        "Addressbook",
        "Sign",
        "Forward",
        "Verify",
        "Decrypt"
    ]

    assert count_interactions_py(sample, concrete_features) == 120
    assert count_interactions(sample, concrete_features) == 120
