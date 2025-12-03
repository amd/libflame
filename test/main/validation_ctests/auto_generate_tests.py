#!/usr/bin/env python3
"""
Generate CMake CTest definitions from YAML test parameters.

Usage:
    python3 auto_generate_tests.py [--all | --api APINAME] <output_dir>

Examples:
    python3 auto_generate_tests.py --api getrf /path/to/output
    python3 auto_generate_tests.py --all /path/to/output
"""

import yaml
import argparse
import sys
from pathlib import Path

def validate_api_definition(api_name, api_def):
    """Validate that API definition has required keys."""
    required_keys = ['param_order', 'description', 'api_group']
    missing_keys = [key for key in required_keys if key not in api_def]
    if missing_keys:
        raise ValueError(f"API '{api_name}' missing required keys: {', '.join(missing_keys)}")
    return True

def validate_path_definition(path, path_index, param_order):
    """Validate that a test path definition has required keys and parameters."""
    if 'name' not in path:
        raise ValueError(f"Path at index {path_index} missing required key: 'name'")
    if 'params' not in path:
        raise ValueError(f"Path '{path.get('name', f'index_{path_index}')}' missing required key: 'params'")

    params = path['params']
    missing_params = [p for p in param_order if p not in params]
    if missing_params:
        raise ValueError(f"Path '{path['name']}' missing required parameters: {', '.join(missing_params)}")
    return True

def load_global_config(script_dir):
    """Load global configuration from auto_generate_label_groups.yaml."""
    label_groups_file = script_dir / "auto_generate_label_groups.yaml"
    config = {
        'label_groups': {},
        'size_thresholds': None,
        'auto_assign_groups': {}
    }

    if label_groups_file.exists():
        try:
            with open(label_groups_file, 'r') as f:
                global_data = yaml.safe_load(f)
                if global_data:
                    config['label_groups'] = global_data.get('label_groups', {})
                    config['size_thresholds'] = global_data.get('size_thresholds')
                    config['auto_assign_groups'] = global_data.get('auto_assign_groups', {})
        except (yaml.YAMLError, IOError) as e:
            print(f"Warning: Failed to load auto_generate_label_groups.yaml: {e}", file=sys.stderr)

    return config

def evaluate_size_label(params, size_thresholds):
    """Evaluate size label based on parameter values and thresholds."""
    if not size_thresholds:
        return None

    param_name = size_thresholds.get('parameter', 'max')
    small_threshold = size_thresholds.get('small')
    medium_threshold = size_thresholds.get('medium')

    if small_threshold is None or medium_threshold is None:
        return None

    # Get the value to check
    if param_name == 'max':
        # Check max of m and n if both exist
        m_val = params.get('m')
        n_val = params.get('n')
        if m_val is not None and n_val is not None:
            value = max(m_val, n_val)
        elif n_val is not None:
            value = n_val
        elif m_val is not None:
            value = m_val
        else:
            return None
    elif param_name in params:
        value = params[param_name]
    else:
        return None

    # Evaluate thresholds
    try:
        value = int(value)
        if value < small_threshold:
            return 'small'
        elif value < medium_threshold:
            return 'medium'
        else:
            return 'large'
    except (ValueError, TypeError):
        return None

def auto_assign_group_labels(size_label, auto_assign_groups):
    """Auto-assign group labels based on size label."""
    group_labels = []
    for group_name, size_labels in auto_assign_groups.items():
        if size_label in size_labels:
            group_labels.append(group_name)
    return group_labels

def expand_label_groups(label_list, label_groups):
    """Expand label groups to their constituent labels. Group label is also included."""
    expanded = []
    for label in label_list:
        if label in label_groups:
            expanded.append(label)
            expanded.extend(label_groups[label])
        else:
            expanded.append(label)
    return expanded

def generate_cmake_for_api(api_name, api_def, output_file, global_config=None):
    """Generate CMake test definitions for a single API."""
    # Validate API definition
    validate_api_definition(api_name, api_def)

    param_order = api_def['param_order']
    code_paths = api_def.get('custom_paths', [])

    # Get global configuration
    if global_config is None:
        global_config = {
            'label_groups': {},
            'size_thresholds': None,
            'auto_assign_groups': {}
        }

    label_groups = global_config['label_groups']
    size_thresholds = global_config['size_thresholds']
    auto_assign_groups = global_config['auto_assign_groups']

    with open(output_file, 'w') as out:
        # File header
        out.write(f"# Auto-generated CTest definitions for {api_name.upper()}\n")
        out.write(f"# DO NOT EDIT - Generated from {api_name}.yaml\n\n")

        # Generate tests
        for idx, path in enumerate(code_paths):
            # Validate path definition
            validate_path_definition(path, idx, param_order)

            path_name = path['name']
            params = path['params']

            # Get precision string (supports both list and string format)
            prec_raw = path.get('precisions', api_def.get('precisions', []))
            if not prec_raw:
                raise ValueError(f"Path '{path_name}' missing precisions (neither in path nor API definition)")
            prec_string = ''.join(prec_raw) if isinstance(prec_raw, list) else prec_raw

            # Build test name and command
            # Match existing naming convention: API name UPPERCASE, precision and path lowercase
            # Examples: GETRF_sd_kernel_2x2, SYEV_sdcz_single_element, LU_FACTORIZATION_sgetrf_2x2
            test_name = f"{api_name.upper()}_{prec_string}_{path_name}"
            cmd_params = ' '.join(str(params[p]) for p in param_order)

            # Build labels from YAML (path-level takes precedence, then API-level, then defaults)
            labels = []

            # Check if path has labels and if they're non-empty
            path_has_labels = 'labels' in path
            path_labels_list = []
            if path_has_labels:
                path_labels_list = path['labels'] if isinstance(path['labels'], list) else [path['labels']]
                path_has_labels = len(path_labels_list) > 0

            # Use API-level labels if path doesn't have labels or has empty labels
            if not path_has_labels and 'labels' in api_def:
                api_labels = api_def['labels'] if isinstance(api_def['labels'], list) else [api_def['labels']]
                labels.extend(expand_label_groups(api_labels, label_groups))

            # Add path-specific labels if present and non-empty (takes precedence over API-level)
            if path_has_labels:
                labels.extend(expand_label_groups(path_labels_list, label_groups))

            # Always add default labels
            # Add both uppercase and lowercase versions for case-insensitive matching
            # This allows users to query with either "GETRF" or "getrf" and both will work
            api_group = api_def['api_group']
            labels.extend([
                api_name.upper(),      # e.g., "GETRF" - allows ctest -L "GETRF"
                api_name.lower(),      # e.g., "getrf" - allows ctest -L "getrf"
                api_group.upper(),     # e.g., "LIN" - allows ctest -L "LIN"
                api_group.lower(),     # e.g., "lin" - allows ctest -L "lin"
                f"precision_{prec_string}"
            ])

            # Add individual precision labels for filtering (e.g., precision_d, precision_s)
            precision_map = {'s': 'precision_s', 'd': 'precision_d', 'c': 'precision_c', 'z': 'precision_z'}
            for prec_char in prec_string:
                if prec_char in precision_map:
                    prec_label = precision_map[prec_char]
                    if prec_label not in labels:
                        labels.append(prec_label)

            # Automatic size-based labeling from parameter thresholds
            size_label = evaluate_size_label(params, size_thresholds)
            if size_label and size_label not in labels:
                labels.append(size_label)

                # Auto-assign group labels based on size label
                auto_group_labels = auto_assign_group_labels(size_label, auto_assign_groups)
                for group_label in auto_group_labels:
                    if group_label not in labels:
                        labels.append(group_label)

            # Add optional arch/size tags based on path name (if not already in labels)
            # Note: avx2/avx512 labels should be set explicitly in YAML, not auto-detected
            # This is a fallback if size_thresholds are not configured
            if not size_label and any(size in path_name for size in ['small', 'medium', 'large']):
                path_size_label = next((s for s in ['small', 'medium', 'large'] if s in path_name), None)
                if path_size_label and path_size_label not in labels:
                    labels.append(path_size_label)

            # Remove duplicates while preserving order
            seen = set()
            labels = [x for x in labels if not (x in seen or seen.add(x))]

            # Write test definition
            out.write(f"add_test(NAME {test_name}\n")
            out.write(f"         COMMAND ${{CTEST_MAIN_COMMAND}} {api_name} {prec_string} {cmd_params}\n")
            out.write(f"         WORKING_DIRECTORY ${{CTEST_WORKING_DIR}})\n")
            out.write(f"set_tests_properties({test_name} PROPERTIES\n")
            out.write(f"    LABELS \"{';'.join(labels)}\"\n")
            out.write(f"    FAIL_REGULAR_EXPRESSION \"FAIL;No test was run\")\n\n")

def main():
    parser = argparse.ArgumentParser(description='Generate CTests from YAML parameters')
    parser.add_argument('--api', help='Generate tests for specific API')
    parser.add_argument('--all', action='store_true', help='Generate tests for all APIs')
    parser.add_argument('output_dir', help='Output directory for generated CMake test files')
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load global configuration (label groups, size thresholds, auto-assign groups)
    global_config = load_global_config(script_dir)

    if args.api:
        yaml_file = script_dir / f"{args.api}.yaml"
        if not yaml_file.exists():
            print(f"Error: {yaml_file} not found", file=sys.stderr)
            return 1

        try:
            with open(yaml_file, 'r') as f:
                data = yaml.safe_load(f)
        except yaml.YAMLError as e:
            print(f"Error: Failed to parse YAML file {yaml_file}: {e}", file=sys.stderr)
            return 1
        except IOError as e:
            print(f"Error: Failed to read file {yaml_file}: {e}", file=sys.stderr)
            return 1

        if data is None:
            print(f"Error: YAML file {yaml_file} is empty or invalid", file=sys.stderr)
            return 1

        api_name = args.api
        if api_name not in data:
            print(f"Error: API '{api_name}' not found in {yaml_file}", file=sys.stderr)
            return 1

        api_def = data[api_name]
        output_file = output_dir / f"{api_name}_tests.cmake"

        try:
            generate_cmake_for_api(api_name, api_def, output_file, global_config)
            print(f"Generated {len(api_def.get('custom_paths', []))} tests for {api_name.upper()}")
        except (ValueError, KeyError) as e:
            print(f"Error: Failed to generate tests for {api_name}: {e}", file=sys.stderr)
            return 1

    elif args.all:
        print("Generating tests for all APIs...")
        errors = []
        for yaml_file in script_dir.glob('*.yaml'):
            # Skip auto_generate_label_groups.yaml as it's not an API definition file
            if yaml_file.stem == 'auto_generate_label_groups':
                continue
            api_name = yaml_file.stem
            try:
                with open(yaml_file, 'r') as f:
                    data = yaml.safe_load(f)
            except yaml.YAMLError as e:
                errors.append(f"  {api_name.upper()}: YAML parse error - {e}")
                continue
            except IOError as e:
                errors.append(f"  {api_name.upper()}: File read error - {e}")
                continue

            if data is None:
                errors.append(f"  {api_name.upper()}: Empty or invalid YAML file")
                continue

            if api_name not in data:
                errors.append(f"  {api_name.upper()}: API not found in YAML")
                continue

            api_def = data[api_name]
            output_file = output_dir / f"{api_name}_tests.cmake"
            try:
                generate_cmake_for_api(api_name, api_def, output_file, global_config)
                print(f"  {api_name.upper()}: {len(api_def.get('custom_paths', []))} tests")
            except (ValueError, KeyError) as e:
                errors.append(f"  {api_name.upper()}: {e}")
                continue

        if errors:
            print("\nErrors encountered:", file=sys.stderr)
            for error in errors:
                print(error, file=sys.stderr)
            return 1
    else:
        parser.print_help()
        return 1

    return 0

if __name__ == '__main__':
    exit(main())
